/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Synthetic five-body benchmark. AU, yr, solar masses; no ephemeris claim. */
#include "rebound_bridge.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double G = REB_BRIDGE_G_AU_YR_MSUN;
static const double masses[5] = {1., 1.e-3, 5.e-5, 8.e-4, 4.e-5};
static const double internal_a[2] = {.02, .025};
struct state {
    struct reb_simulation* main;
    struct reb_simulation* sub[2];
    struct reb_bridge* bridge;
};

static void require(int condition, const char* message) {
    if (!condition) {
        fprintf(stderr, "%s: %s\n", message, reb_bridge_last_error_message());
        exit(EXIT_FAILURE);
    }
}

static struct reb_simulation* sim_new(double epsilon, double dt) {
    struct reb_simulation* r = reb_simulation_create();
    require(r != NULL, "allocation failed");
    r->G = G;
    require(reb_bridge_set_integrator_ias15(r) == 0, "IAS15 setup");
    r->ri_ias15.epsilon = epsilon;
    r->exact_finish_time = 1;
    r->dt = dt;
    return r;
}

static void get_com(const struct reb_particle* p, int first, double* R, double* V) {
    const double m = p[first].m + p[first+1].m;
    R[0]=(p[first].m*p[first].x+p[first+1].m*p[first+1].x)/m;
    R[1]=(p[first].m*p[first].y+p[first+1].m*p[first+1].y)/m;
    R[2]=(p[first].m*p[first].z+p[first+1].m*p[first+1].z)/m;
    V[0]=(p[first].m*p[first].vx+p[first+1].m*p[first+1].vx)/m;
    V[1]=(p[first].m*p[first].vy+p[first+1].m*p[first+1].vy)/m;
    V[2]=(p[first].m*p[first].vz+p[first+1].m*p[first+1].vz)/m;
}

static struct reb_particle proxy(const struct reb_particle* p, int first) {
    double R[3],V[3]; get_com(p,first,R,V);
    struct reb_particle q={0};
    q.m=p[first].m+p[first+1].m;
    q.x=R[0];q.y=R[1];q.z=R[2];q.vx=V[0];q.vy=V[1];q.vz=V[2];
    return q;
}

static void initial_particles(struct reb_particle* p) {
    memset(p,0,5*sizeof(*p));p[0].m=masses[0];
    const double outer_a[2]={1.,1.45};
    const double phase[2]={.3,1.8};
    const double ecc[2]={.05,.04}, inc[2]={.1,-.12};
    for(int s=0;s<2;s++) {
        const int first=1+2*s;
        struct reb_simulation* r=sim_new(1.e-13,1.e-5);
        reb_simulation_add_fmt(r,"m",masses[first]);
        reb_simulation_add_fmt(r,"m a e inc f",masses[first+1],internal_a[s],ecc[s],inc[s],phase[s]);
        reb_simulation_move_to_com(r);
        for(int i=0;i<2;i++) {
            p[first+i]=r->particles[i];
            p[first+i].x+=outer_a[s];
            p[first+i].vy+=sqrt(G*(masses[0]+masses[first]+masses[first+1])/outer_a[s]);
        }
        reb_simulation_free(r);
    }
    struct reb_simulation* full=sim_new(1.e-13,1.e-5);
    for(int i=0;i<5;i++)reb_simulation_add(full,p[i]);
    reb_simulation_move_to_com(full);
    for(int i=0;i<5;i++) {p[i]=full->particles[i];p[i].sim=NULL;}
    reb_simulation_free(full);
}

static void point_acc(const struct reb_particle* target, const struct reb_particle* source, double* a) {
    double dr[3]={source->x-target->x,source->y-target->y,source->z-target->z};
    double r2=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
    double f=G*source->m/(r2*sqrt(r2));
    for(int k=0;k<3;k++)a[k]+=f*dr[k];
}

/* Independently evaluate -grad(U_s0 + U_0t - U_00) in physical coordinates.
 * Differentiating a proxy position distributes its reaction to every member. */
static void approximate_cross_acc(const struct reb_particle* p, double a[5][3]) {
    struct reb_particle ca=proxy(p,1),cb=proxy(p,3);
    double at_com_a[3]={0},at_com_b[3]={0},mono_a[3]={0},mono_b[3]={0};
    for(int j=3;j<5;j++)point_acc(&ca,&p[j],at_com_a);
    for(int i=1;i<3;i++)point_acc(&cb,&p[i],at_com_b);
    point_acc(&ca,&cb,mono_a);point_acc(&cb,&ca,mono_b);
    for(int i=1;i<3;i++) {
        point_acc(&p[i],&cb,a[i]);
        for(int k=0;k<3;k++)a[i][k]+=at_com_a[k]-mono_a[k];
    }
    for(int j=3;j<5;j++) {
        point_acc(&p[j],&ca,a[j]);
        for(int k=0;k<3;k++)a[j][k]+=at_com_b[k]-mono_b[k];
    }
}

static void exact_cross_acc(const struct reb_particle* p, double a[5][3]) {
    for(int i=1;i<3;i++)for(int j=3;j<5;j++) {
        point_acc(&p[i],&p[j],a[i]);point_acc(&p[j],&p[i],a[j]);
    }
}

/* IAS15 reference of the approximate Hamiltonian: full Newtonian force plus
 * an independently evaluated approximate-minus-exact inter-subsystem force. */
static void approximate_correction(struct reb_simulation* r) {
    double a[5][3]={{0}},b[5][3]={{0}};
    approximate_cross_acc(r->particles,a);exact_cross_acc(r->particles,b);
    for(int i=1;i<5;i++) {
        r->particles[i].ax+=a[i][0]-b[i][0];
        r->particles[i].ay+=a[i][1]-b[i][1];
        r->particles[i].az+=a[i][2]-b[i][2];
    }
}

static struct state bridge_new(const struct reb_particle* p, double h, double epsilon,
                               int resolved, int lf42) {
    struct state s={0};s.main=sim_new(epsilon,h/8.);
    reb_simulation_add(s.main,p[0]);
    for(int j=0;j<2;j++)reb_simulation_add(s.main,proxy(p,1+2*j));
    s.bridge=reb_bridge_create(s.main,h);require(s.bridge!=NULL,"create bridge");
    reb_bridge_set_owns_main(s.bridge,1);
    require(reb_bridge_set_subsystem_coupling(s.bridge,resolved?REB_BRIDGE_SUBSYSTEM_COUPLING_RESOLVED:REB_BRIDGE_SUBSYSTEM_COUPLING_MONOPOLE)==0,"coupling");
    require(reb_bridge_set_method(s.bridge,lf42?REB_BRIDGE_METHOD_EOS_LF4_2:REB_BRIDGE_METHOD_ORDER2)==0,"method");
    for(int j=0;j<2;j++) {
        s.sub[j]=sim_new(epsilon,h/8.);
        struct reb_particle c=proxy(p,1+2*j);
        for(int i=0;i<2;i++) {
            struct reb_particle q=p[1+2*j+i];
            q.x-=c.x;q.y-=c.y;q.z-=c.z;q.vx-=c.vx;q.vy-=c.vy;q.vz-=c.vz;
            reb_simulation_add(s.sub[j],q);
        }
        /* Roundoff-only recentering; main proxy retains the original COM. */
        reb_simulation_move_to_com(s.sub[j]);
        require(reb_bridge_add_subsystem(s.bridge,j+1,s.sub[j],h/8.,1)==0,"register");
    }
    return s;
}

static void reconstruct(const struct state* s, struct reb_particle* p) {
    p[0]=s->main->particles[0];
    for(int j=0;j<2;j++)for(int i=0;i<2;i++) {
        p[1+2*j+i]=s->sub[j]->particles[i];
        struct reb_particle* q=&p[1+2*j+i];
        const struct reb_particle* c=&s->main->particles[j+1];
        q->x+=c->x;q->y+=c->y;q->z+=c->z;
        q->vx+=c->vx;q->vy+=c->vy;q->vz+=c->vz;
    }
}

static void output(FILE* f,double days,const struct reb_particle* p) {
    for(int i=0;i<5;i++)fprintf(f,"%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
        days,i,p[i].m,p[i].x,p[i].y,p[i].z,p[i].vx,p[i].vy,p[i].vz);
}

static int kick_check(const char* path) {
    struct reb_particle initial[5];initial_particles(initial);
    FILE* f=fopen(path,"w");require(f!=NULL,"open kick output");
    fprintf(f,"coupling,body,component,actual_B,expected_B,abs_error\n");
    double maxerr=0.;
    for(int resolved=0;resolved<2;resolved++) {
        struct state s=bridge_new(initial,.01,1.e-13,resolved,0);
        struct reb_particle before[5],after[5];reconstruct(&s,before);
        double expected[5][3]={{0}},drift[5][3]={{0}};
        /* Whole physical Hamiltonian force, independently constructed. */
        for(int i=0;i<5;i++)for(int j=0;j<5;j++)if(i!=j)point_acc(&before[i],&before[j],expected[i]);
        if(!resolved) {
            double a[5][3]={{0}},b[5][3]={{0}};
            approximate_cross_acc(before,a);exact_cross_acc(before,b);
            for(int i=0;i<5;i++)for(int k=0;k<3;k++)expected[i][k]+=a[i][k]-b[i][k];
        }
        double main_a[3][3]={{0}};
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(i!=j)point_acc(&s.main->particles[i],&s.main->particles[j],main_a[i]);
        for(int k=0;k<3;k++)drift[0][k]=main_a[0][k];
        for(int j=0;j<2;j++)for(int i=0;i<2;i++) {
            const int idx=1+2*j+i;
            point_acc(&before[idx],&before[1+2*j+1-i],drift[idx]);
            for(int k=0;k<3;k++)drift[idx][k]+=main_a[j+1][k];
        }
        const double dt=.125;
        require(reb_bridge_apply_cross_kick(s.bridge,dt)==0,"kick");reconstruct(&s,after);
        for(int i=0;i<5;i++) {
            double v0[3]={before[i].vx,before[i].vy,before[i].vz};
            double v1[3]={after[i].vx,after[i].vy,after[i].vz};
            for(int k=0;k<3;k++) {
                double actual=(v1[k]-v0[k])/dt, want=expected[i][k]-drift[i][k];
                double err=fabs(actual-want);if(err>maxerr)maxerr=err;
                fprintf(f,"%s,%d,%d,%.17g,%.17g,%.17g\n",resolved?"resolved":"monopole",i,k,actual,want,err);
            }
        }
        reb_bridge_free(s.bridge);
    }
    fclose(f);printf("kick_max_abs_error=%.17g\n",maxerr);
    require(maxerr<1.e-10,"kick differs from independent Hamiltonian force");return 0;
}

int main(int argc,char** argv) {
    if(argc==3 && strcmp(argv[1],"--kick-check")==0)return kick_check(argv[2]);
    if(argc!=9) {
        fprintf(stderr,"Usage: %s bridge|reference monopole|resolved kdk|lf4_2 h_days epsilon duration_days sample_days output.csv\n",argv[0]);
        return EXIT_FAILURE;
    }
    int is_bridge=strcmp(argv[1],"bridge")==0,resolved=strcmp(argv[2],"resolved")==0;
    int lf42=strcmp(argv[3],"lf4_2")==0;
    require(is_bridge || strcmp(argv[1],"reference")==0,"invalid mode");
    require(resolved || strcmp(argv[2],"monopole")==0,"invalid coupling");
    require(lf42 || strcmp(argv[3],"kdk")==0,"invalid method");
    double hdays=atof(argv[4]),epsilon=atof(argv[5]),duration=atof(argv[6]),sample=atof(argv[7]);
    require(hdays>0 && epsilon>0 && duration>0 && sample>0,"invalid positive argument");
    int n=(int)llround(duration/sample),steps=(int)llround(sample/hdays);
    require(fabs(n*sample-duration)<1.e-10 && steps>0 && fabs(steps*hdays-sample)<1.e-10,"sample alignment");
    struct reb_particle p[5];initial_particles(p);
    FILE* f=fopen(argv[8],"w");require(f!=NULL,"open output");
    fprintf(f,"t_days,body,m,x,y,z,vx,vy,vz\n");output(f,0,p);
    struct state s={0};struct reb_simulation* reference=NULL;
    if(is_bridge)s=bridge_new(p,hdays/365.25,epsilon,resolved,lf42);
    else {
        reference=sim_new(epsilon,hdays/365.25/8.);
        for(int i=0;i<5;i++)reb_simulation_add(reference,p[i]);
        if(!resolved)reference->additional_forces=approximate_correction;
    }
    for(int i=1;i<=n;i++) {
        const double target=i*sample/365.25;
        if(is_bridge) {
            /* Public entry point; all output intervals contain integral outer steps. */
            require(reb_bridge_advance(s.bridge,sample/365.25)==0,"advance");
            reconstruct(&s,p);
            require(fabs(s.main->t-target)<1.e-10,"sample epoch");
        }else {
            require(reb_simulation_integrate(reference,target)==REB_STATUS_SUCCESS,"reference advance");
            for(int j=0;j<5;j++)p[j]=reference->particles[j];
        }
        output(f,i*sample,p);
    }
    fclose(f);
    if(is_bridge)reb_bridge_free(s.bridge);else reb_simulation_free(reference);
    return EXIT_SUCCESS;
}
