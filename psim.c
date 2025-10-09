#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL2/SDL.h>
#include <sys/resource.h>
#include <time.h>

const char control_fname[] = "/home/jtmcg/projects/psim/render_color.txt";

// Window geometry
#define WIDTH 5500.0
#define RENDER_WIDTH 3680.0
#define HEIGHT 1000.0
#define WALL_RADIUS 4500.0
#define NECK_HEIGHT 250.0

// Simulation parameters
#define NUM_PARTICLES 4500 //number of particles in the simulation
#define RADIUS 5.5 //radius of each particle
#define PPS 400 // Target average speed in channel in pixels per second for viewing
#define DT 1.0 // initial timestep for simulation
#define REPOS 1.0001 // Repositioning constant. Puts particles slightly further away on bounce to prevent them getting stuck

// Blower power, expressed as a fraction of initial temperature
#define BLOW_POWER 1.0
 
// Enable this to apply a periodic acceleration to every particle instead of a "fan" model
//     where acceleration is applied as the particles cross the end plane
//#define MAGIC_FORCE
//#define ADIABATIC_CHANNEL

// Surface roughness, scales the random shifts applied to wall bounce angle
//     Effectively generates entropy
#define ROUGHNESS_ANG 0.0

// Number of distinct sections for pressure and speed computations
#define SECTIONS 18
#define PFRAMES 500
#define KICK_PERIOD 20

// This is a crude approximation for intermolecular forces
#define FORCE 0.0
#define FPL_MULT 0.8

//#define INIT_HALF

// Define this to draw walls. Takes a while on startup to compute the draw circle
//#define DRAW_WALLS
const int wall_color[3] = {128, 128, 140};


//structure for a single particle
typedef struct {
    float x; //x position of particle
    float y; //y position of particle
    float vx; //x velocity of particle
    float vy; //y velocity of particle
} particle;

// Returns a random float between 0 and 1
inline float randf(){ return ((float)rand() / RAND_MAX); }

static float dt = DT;

// Function to read information from control file
char get_render_color(float *blow){
    FILE* fp = fopen(control_fname, "r");
    if( fp == NULL ){
        return 0;
    }

    int char0 = fgetc(fp);
    int charr = 0;
    float res = 0.0;
    float dmul = 1000.0;
    uint8_t valid = 0;

    while (char0!=EOF && char0!='\n'){
        if(char0>=48 && char0<=57){
            dmul = dmul/10.0;
            res += (float)(char0-48)*dmul;
        }
        else if(char0 == '.'){
            res = res/dmul;
            dmul = 1.0;
            valid = 1;
        }else if(char0!=',' && char0!='\n' && char0!=' ' && char0!='\t'){
            charr = char0;
        }
        char0 = fgetc(fp);
    }

    if(valid){
        *blow = res;
    }else{
        *blow = BLOW_POWER;
    }

    fclose(fp);
    return (char)charr;
}

// Function to randomize the input bounce angle according to the roughness angle
void rand_bounce(float *nx, float *ny){
    float rang = (randf()-0.5)*2*ROUGHNESS_ANG;
    float c = cosf(rang);
    float s = sinf(rang);
    float ox = *nx;
    float oy = *ny;
    *nx = ox*c - oy*s;
    *ny = ox*s + oy*c;
}

// Householder reflection function for particle velocities
void householder(particle *p, float nx, float ny, float bleed){
    float s = (1.0+bleed)*(nx*p->vx + ny*p->vy);
    p->vx -= nx*s;
    p->vy -= ny*s;
}

// Geometry definitions for walls
#define BIGCIRC_X (WIDTH/2)
#define UNREND ((WIDTH-RENDER_WIDTH)/2.0)
#define BIGCIRC_YH (HEIGHT/2+NECK_HEIGHT/2+WALL_RADIUS)
#define BIGCIRC_YL (HEIGHT/2-NECK_HEIGHT/2-WALL_RADIUS)
#define PUSH_RAD ((WALL_RADIUS+RADIUS)*REPOS)

static float bounce_accum[SECTIONS];
static float speeds[SECTIONS];
static float x_spds[SECTIONS];
static float pressure[SECTIONS];
static float pressure_igl[SECTIONS];
static float density[SECTIONS];
static float temp[SECTIONS];
static float temperature[SECTIONS];
static int spd_cnt[SECTIONS];
static float p_areas[SECTIONS];
static float d_areas[SECTIONS];
static double speed_sum = 0;
static int64_t speed_cnt = 0;

static float mean_fpl;

#define BIGCIRC_XR (BIGCIRC_X-UNREND)
void get_areas(){
    float lb, rb, cx, cy, ar, nar;
    float dx = sqrt(WALL_RADIUS*WALL_RADIUS - BIGCIRC_YL*BIGCIRC_YL);
    float theta0 = asinf(dx/WALL_RADIUS);
    float theta1 = asinf(-BIGCIRC_YL/WALL_RADIUS);
    float lc_bound = BIGCIRC_XR - dx;
    float rc_bound = BIGCIRC_XR + dx;
    for(int i=0; i<SECTIONS; i++){
        // Get left/right x bounds
        lb = (RENDER_WIDTH/SECTIONS)*i;
        rb = lb + (RENDER_WIDTH/SECTIONS);
        // Left of the big circle
        if(lb<lc_bound && rb<=lc_bound){
            ar = rb-lb;
            nar = HEIGHT;
        // Bounds contain the left edge of the big circle
        }else if(rb>lc_bound && lb<lc_bound){
            float theta = asinf((BIGCIRC_XR-rb)/WALL_RADIUS);
            ar = (lc_bound-lb) + WALL_RADIUS*(theta0 - theta);
            nar = (HEIGHT + HEIGHT - 2*(BIGCIRC_YL+WALL_RADIUS*cosf(theta)))/2;
        // Bounds within big circle
        }else if(lb>=lc_bound && rb<=rc_bound){
            float thetaA = asinf((BIGCIRC_XR-lb)/WALL_RADIUS);
            float thetaB = asinf((BIGCIRC_XR-rb)/WALL_RADIUS);
            ar = WALL_RADIUS * (thetaA - thetaB);
            nar = (HEIGHT-2*(BIGCIRC_YL+WALL_RADIUS*cosf(thetaA)) + HEIGHT-2*(BIGCIRC_YL+WALL_RADIUS*cosf(thetaB)))/2;
        // Bounds contain the right edge of the big circle
        }else if(lb<rc_bound && rb>rc_bound){
            float theta = asinf((BIGCIRC_XR-lb)/WALL_RADIUS);
            ar = WALL_RADIUS*(theta + theta0) + (rb-rc_bound);
            nar = (HEIGHT + HEIGHT - 2*(BIGCIRC_YL+WALL_RADIUS*cosf(theta)))/2;
        // Otherwise assume right of big circle
        }else{
            ar = rb-lb;
            nar = HEIGHT;
        }
        
        p_areas[i] = 2*ar;
        d_areas[i] = nar;
        printf("SECTION %d: WALL AREA %.1f, NECK AREA %.1f\n", i, 2*ar, nar);
    }
}

static float init_avg_speed = 0;
static float avg_speed = 0;
static float bpwr = BLOW_POWER, bpwr_command = BLOW_POWER;

inline uint8_t top_bot_bounce(particle *p, float energy_bleed, const float wall_offset){
    float py = p->y;
    // Check Y rectangle bounds
    if(py <= (wall_offset + RADIUS)){
        householder(p, 0, -1, energy_bleed);
        p->y = RADIUS*REPOS + wall_offset;
        return 1;
    }
    if(py >= (HEIGHT-RADIUS-wall_offset)){
        householder(p, 0, 1, energy_bleed);
        p->y = HEIGHT - RADIUS*REPOS - wall_offset;
        return 1;
    }
    return 0;
}

inline uint8_t constriction(particle *p, float energy_bleed, const uint8_t open_close, const float throat_width, const float length, const float throat_x){
    // Distance from wall to throat opening
    const float b = (HEIGHT-throat_width)/2;
    const float m = b / length;
    const float m2 = m*m;
    const float nx = sqrtf(m2 / (1+m2));
    const float ny = sqrtf(1.0 / (1+m2));
    if(open_close){
        // Return no change if outside of x bounds
        //if((p->x < (throat_x-length)) || (p->x >= throat_x))  return 0;
        // Compute slope and y intercept of bounding lines
        const float b_up = m*(length-throat_x) + RADIUS;
        float yc = m*(p->x) + b_up;
        // Check upper boundary
        if(p->y < yc){
            householder(p, -nx, ny, energy_bleed);
            p->y = yc+(REPOS-1.0);
            return 1;
        }
        // Repeat for lower boundary
        const float b_dwn = (-m)*(length-throat_x) + HEIGHT - RADIUS;
        yc = (-m)*(p->x) + b_dwn;
        if(p->y > yc){
            householder(p, -nx, -ny, energy_bleed);
            p->y = yc-(REPOS-1.0);
            return 1;
        }
    }else{
        // Return no change if outside of x bounds
        //if((p->x < throat_x) || (p->x >= (throat_x+length)))  return 0;
        // Compute slope and y intercept of bounding lines
        const float b_up = m*throat_x + RADIUS + b;
        float yc = (-m)*(p->x) + b_up;
        // Check upper boundary
        if(p->y < yc){
            householder(p, nx, ny, energy_bleed);
            p->y = yc+(REPOS-1.0);
            return 1;
        }
        // Repeat for lower boundary
        const float b_dwn = (-m)*throat_x + b + throat_width - RADIUS;
        yc = m*(p->x) + b_dwn;
        if(p->y > yc){
            householder(p, nx, -ny, energy_bleed);
            p->y = yc-(REPOS-1.0);
            return 1;
        }
    }
    return 0;
}

inline int get_section(particle *p){
    float s = floorf((p->x - UNREND)/(RENDER_WIDTH/SECTIONS));
    int si = (int)s;
    if(si>(SECTIONS-1))  si = -1;
    return si;
}

#define THROAT_WIDTH (HEIGHT*0.3)
uint8_t venturi_bounce2(particle *p, float energy_bleed){
    // Implement left boundary wrap with magic force optional
    if(p->x <= 0){
        p->x = (WIDTH-REPOS+1);
        #ifndef MAGIC_FORCE
            p->vx -= bpwr_command;
        #endif
    // Implement right boundary crossing
    }else if(p->x >= WIDTH)  p->x = (REPOS-1);
    
    int s = get_section(p);
    float vx = p->vx, vy = p->vy;
    speed_sum += vx;
    speed_cnt++;

    if(s > -1){
        x_spds[s] += vx;
        temperature[s] += (vx*vx + vy*vy);
        spd_cnt[s] ++;
    }
    // Compute the X offset
    const float x_offset = (WIDTH - RENDER_WIDTH)*0.5;
    float x = p->x - x_offset;
    // Wide wall bounce
    if( x < RENDER_WIDTH*0.15 ){
        return top_bot_bounce(p, energy_bleed, 0.0);
    // Divergent section
    }else if( x < RENDER_WIDTH*0.5 ){
        return constriction(p, energy_bleed, 1, THROAT_WIDTH, RENDER_WIDTH*(0.5-0.15), RENDER_WIDTH*0.5+x_offset);
    }else if( x < RENDER_WIDTH*0.65 ){
        return top_bot_bounce(p, energy_bleed, (HEIGHT-THROAT_WIDTH)*0.5);
    }else if( x < RENDER_WIDTH*0.85 ){
        return constriction(p, energy_bleed, 0, THROAT_WIDTH, RENDER_WIDTH*(0.85-0.65), RENDER_WIDTH*0.65+x_offset);
    }else{
        return top_bot_bounce(p, energy_bleed, 0.0);
    }
}

uint8_t venturi_bounce(particle *p, float energy_bleed){
    float px = p->x;
    // Check X rectangle bounds
    if(px <= 0){
        p->x = (WIDTH-REPOS+1);
        #ifndef MAGIC_FORCE
            p->vx -= bpwr_command;
        #endif
        //p->vy *= BLOW_POWER/2;
        return 1;
    }
    if(px >= WIDTH){
        p->x = (REPOS-1);
        return 1;
    }

    px = p->x;
    int pind = (int)((px - UNREND)/(RENDER_WIDTH/SECTIONS));
    uint8_t pvalid = (px < RENDER_WIDTH+UNREND) && (px >= UNREND);
    float vx = p->vx, vy = p->vy;
    speed_sum += vx;
    speed_cnt++;

    if(pvalid){
        x_spds[pind] += vx;
        temperature[pind] += (vx*vx + vy*vy);
        spd_cnt[pind]++;
    }

    float s, nx, ny;
    float py = p->y;
    // Check Y rectangle bounds
    if(py <= RADIUS){
        float ay = p->vy;
        if(ay<0){ ay= -ay; }
        if(pvalid){
            bounce_accum[pind] += ay;
        }
        nx=0;
        ny=-1;
        rand_bounce(&nx, &ny);
        householder(p, nx, ny, energy_bleed);
        p->y = RADIUS*REPOS;
        return 1;
    }
    if(py >= (HEIGHT-RADIUS)){
        float ay = p->vy;
        if(ay<0){ ay= -ay; }
        if(pvalid){
            bounce_accum[pind] += ay;
        }
        nx=0;
        ny=1;
        rand_bounce(&nx, &ny);
        householder(p, nx, ny, energy_bleed);
        p->y = HEIGHT-RADIUS*REPOS;
        return 1;
    }

    // Compute distance from upper part of constriction bound
    float dx = p->x - BIGCIRC_X;
    float dy = p->y - BIGCIRC_YH;
    float dx2 = dx*dx;
    float dist = sqrt(dx2 + dy*dy);
    float mom;
    if( dist<(WALL_RADIUS+RADIUS) ){
        // Perform householder reflection of velocity
        nx = dx/dist;
        ny = dy/dist;
        mom = nx*p->vx + ny*p->vy; // Scalar projection of the particle velocity onto the normal vector
        if(mom<0){ mom = -mom; }
        bounce_accum[pind] += mom;
        p->x = BIGCIRC_X + PUSH_RAD*nx;
        p->y = BIGCIRC_YH + PUSH_RAD*ny;
        rand_bounce(&nx, &ny);
        #ifdef ADIABATIC_CHANNEL
        householder(p, nx, ny, 1.0);
        #else
        householder(p, nx, ny, energy_bleed);
        #endif
        return 1;
    }

    dy = p->y - BIGCIRC_YL;
    dist = sqrt(dx2 + dy*dy);
    if( dist<(WALL_RADIUS+RADIUS) ){
        // Perform householder reflection of velocity
        nx = dx/dist;
        ny = dy/dist;
        mom = nx*p->vx + ny*p->vy;
        if(mom<0){ mom = -mom; }
        bounce_accum[pind] += mom;
        p->x = BIGCIRC_X + PUSH_RAD*nx;
        p->y = BIGCIRC_YL + PUSH_RAD*ny;
        rand_bounce(&nx, &ny);
        #ifdef ADIABATIC_CHANNEL
        householder(p, nx, ny, 1.0);
        #else
        householder(p, nx, ny, energy_bleed);
        #endif
        return 1;
    }
    return 0;
}

// Function to check if particles overlap with themselves or bounds
inline uint8_t check_overlap(particle* particles, float x, float y, int occ){
    particle p;
    p.x=x; p.y=y; p.vx=0; p.vy=0;
    #ifdef INIT_HALF
    if(x<WIDTH*0.65 || x>WIDTH*0.9){ return 1; }
    #endif
    if(venturi_bounce2(&p, 1.0)){ return 1; }
    for(int i=0; i<occ; i++){
        float dx = particles[i].x - x;
        float dy = particles[i].y - y;
        float dist = sqrt(dx*dx + dy*dy);
        if(dist<=2*RADIUS){ return 1; }
    }
    return 0;
}

// Function to initialize particles with random positions and velocities
void init_particles(particle* particles) {
    int i=0;
    while(i<NUM_PARTICLES){
        float cx = randf()*WIDTH;
        float cy = randf()*HEIGHT;
        if(!check_overlap(particles, cx, cy, i)){
            particles[i].x = cx;
            particles[i].y = cy;
            particles[i].vy = (randf()*2 - 1);
            #ifndef MAGIC_FORCE
                particles[i].vx = (randf()*2 - 1);
            #else
                particles[i].vx = ((randf()*2 - 1) - BLOW_POWER);
            #endif
            i++;
        }
    }
}

// Declare and zero the initial energy
static double initial_energy=0.0;
static double slowdown = 1.0;

// Function to simulate a timestep for each particle
float simulate_timestep(particle* particles, float ideal_dist) {
    float spd2;
    // Compute the initial energy if it is the first go around
    if(initial_energy==0.0){
        for (int i = 0; i < NUM_PARTICLES; i++) {
            spd2 = particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy;
            initial_energy += spd2;
        }
    }

    // Scan through all particle pair combos
    for (int i = 0; i < NUM_PARTICLES - 1; i++) {
        float ix = particles[i].x;
        float iy = particles[i].y;
        float ivx = particles[i].vx;
        float ivy = particles[i].vy;
        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            float jx = particles[j].x;
            float jy = particles[j].y;
            float jvx = particles[j].vx;
            float jvy = particles[j].vy;
            // Compute distance with distance formula
            float dx = ix - jx;
            float dy = iy - jy;
            float dist2 = dx*dx + dy*dy;
            float dist = sqrt(dist2);
            // Compute normal vector components
            float nx = dx/dist;
            float ny = dy/dist;
            if (dist <= 2*RADIUS) {
                // Compute scalar impulse value
                float k = nx*(jvx - ivx) + ny*(jvy - ivy);
                // Add the deltaV to each velocity
                ivx += k*nx;
                ivy += k*ny;
                jvx -= k*nx;
                jvy -= k*ny;
                ix = jx + nx*(2*RADIUS*REPOS);
                iy = jy + ny*(2*RADIUS*REPOS);

                particles[i].vx = ivx;
                particles[i].vy = ivy;
                particles[j].vx = jvx;
                particles[j].vy = jvy;
                particles[i].x = ix;
                particles[i].y = iy;
            }else if(FORCE != 0.0){
                if(dist<mean_fpl){
                    float force = mean_fpl/dist*(FORCE/2);
                    ivx += force*nx;
                    ivy += force*ny;
                    jvx -= force*nx;
                    jvy -= force*ny;

                    particles[i].vx = ivx;
                    particles[i].vy = ivy;
                    particles[j].vx = jvx;
                    particles[j].vy = jvy;
                }
                
            }
        }
    }

    double current_energy = 0;    
    // Update positions of particles
    for (int i = 0; i < NUM_PARTICLES; i++) {
        float x = particles[i].x;
        float y = particles[i].y;
        float vx = particles[i].vx;
        float vy = particles[i].vy;
        x += vx * dt;
        y += vy * dt;
        particles[i].x = x;
        particles[i].y = y;
        spd2 = vx*vx + vy*vy;
        current_energy += spd2;
    }
    
    // Compute the amount of slowdown to prevent excessive energy gain
    if(current_energy > 0){
        slowdown = initial_energy/current_energy;
    }

    // Check for collisions with bounding box
    for (int i = 0; i < NUM_PARTICLES; i++) {
        venturi_bounce2(particles+i, slowdown);
    }

    // Return the current amount of slowdown to maintain energy
    return (float)slowdown;
}

// Define a table to hold pixel offsets for rendering a filled in circle
static SDL_Point circle_lut[(int)((RADIUS+2)*(RADIUS+2)*4)];
static int circle_pixcnt=0;

static SDL_Point wall_high[(int)(WALL_RADIUS*WALL_RADIUS*M_PI/2.0)];
static int wall_cnt_h=0;
static SDL_Point wall_low[(int)(WALL_RADIUS*WALL_RADIUS*M_PI/2.0)];
static int wall_cnt_l=0;

// Function to generate tables of pixels for filled in circles
void generate_circle(SDL_Point *lut, int *pixcnt, int x_off, int y_off, int rad){
    
    uint8_t use_bounds = (x_off!=0 || y_off!=0);

    uint8_t yhigh = y_off>HEIGHT/2;

    if(use_bounds){
        for(int x=0; x<RENDER_WIDTH; x++){
            float xo = x - x_off;
            float und = rad*rad - xo*xo;
            if(und>0){
                float y;
                if(yhigh){
                    y = y_off - sqrtf(und);
                    if(y<HEIGHT){
                        for(int i=(int)y; i<HEIGHT; i++){
                            lut[*pixcnt].x = x;
                            lut[*pixcnt].y = i;
                            (*pixcnt)++;
                        }
                    }
                }else{
                    y = y_off + sqrtf(und);
                    if(y>0){
                        for(int i=0; i<(int)y; i++){
                            lut[*pixcnt].x = x;
                            lut[*pixcnt].y = i;
                            (*pixcnt)++;
                        }
                    }
                }
            }
        }
        return;
    }

    float ang=0.0;
    float xf, yf;
    int x, y;
    float ang_inc = 0.3 / (float)rad;
    int radc = rad*2;
    while(ang<(2*M_PI)){
        xf = cosf(ang);
        yf = sinf(ang);
        ang += ang_inc;
        //printf("Angle is %.3f of %.3f, pixcnt %d\n", ang, 2*M_PI, *pixcnt);
        for(float rf=0; rf<=rad; rf+=0.35){
            x = roundf(xf*rf) + x_off;
            y = roundf(yf*rf) + y_off;
            uint8_t done=0;
            int i = *pixcnt - radc;
            if(i<0){ i=0; }
            while(i<(*pixcnt) && !done){
                done |= ((lut[i].x==x) && (lut[i].y==y));
                i++;
            }
            if(!done){
                lut[(*pixcnt)].x = x;
                lut[(*pixcnt)].y = y;
                (*pixcnt) ++;
            }
        }
    }
}

const uint8_t color_lut[256][3] = {
	{255,0,0},
	{253,1,0},
	{251,3,0},
	{249,5,0},
	{247,7,0},
	{245,9,0},
	{243,11,0},
	{241,13,0},
	{239,15,0},
	{237,17,0},
	{235,19,0},
	{233,21,0},
	{231,23,0},
	{229,25,0},
	{227,27,0},
	{225,29,0},
	{223,31,0},
	{221,33,0},
	{219,35,0},
	{217,37,0},
	{215,39,0},
	{213,41,0},
	{211,43,0},
	{209,45,0},
	{207,47,0},
	{205,49,0},
	{203,51,0},
	{201,53,0},
	{199,55,0},
	{197,57,0},
	{195,59,0},
	{193,61,0},
	{191,63,0},
	{189,65,0},
	{187,67,0},
	{185,69,0},
	{183,71,0},
	{181,73,0},
	{179,75,0},
	{177,77,0},
	{175,79,0},
	{173,81,0},
	{171,83,0},
	{169,85,0},
	{167,87,0},
	{165,89,0},
	{163,91,0},
	{161,93,0},
	{159,95,0},
	{157,97,0},
	{155,99,0},
	{153,101,0},
	{151,103,0},
	{149,105,0},
	{147,107,0},
	{145,109,0},
	{143,111,0},
	{141,113,0},
	{139,115,0},
	{137,117,0},
	{135,119,0},
	{133,121,0},
	{131,123,0},
	{129,125,0},
	{127,127,0},
	{125,129,0},
	{123,131,0},
	{121,133,0},
	{119,135,0},
	{117,137,0},
	{115,139,0},
	{113,141,0},
	{111,143,0},
	{109,145,0},
	{107,147,0},
	{105,149,0},
	{103,151,0},
	{101,153,0},
	{99,155,0},
	{97,157,0},
	{95,159,0},
	{93,161,0},
	{91,163,0},
	{89,165,0},
	{87,167,0},
	{85,169,0},
	{83,171,0},
	{81,173,0},
	{79,175,0},
	{77,177,0},
	{75,179,0},
	{73,181,0},
	{71,183,0},
	{69,185,0},
	{67,187,0},
	{65,189,0},
	{63,191,0},
	{61,193,0},
	{59,195,0},
	{57,197,0},
	{55,199,0},
	{53,201,0},
	{51,203,0},
	{49,205,0},
	{47,207,0},
	{45,209,0},
	{43,211,0},
	{41,213,0},
	{39,215,0},
	{37,217,0},
	{35,219,0},
	{33,221,0},
	{31,223,0},
	{29,225,0},
	{27,227,0},
	{25,229,0},
	{23,231,0},
	{21,233,0},
	{19,235,0},
	{17,237,0},
	{15,239,0},
	{13,241,0},
	{11,243,0},
	{9,245,0},
	{7,247,0},
	{5,249,0},
	{3,251,0},
	{1,253,0},
	{0,255,0},
	{0,253,1},
	{0,251,3},
	{0,249,5},
	{0,247,7},
	{0,245,9},
	{0,243,11},
	{0,241,13},
	{0,239,15},
	{0,237,17},
	{0,235,19},
	{0,233,21},
	{0,231,23},
	{0,229,25},
	{0,227,27},
	{0,225,29},
	{0,223,31},
	{0,221,33},
	{0,219,35},
	{0,217,37},
	{0,215,39},
	{0,213,41},
	{0,211,43},
	{0,209,45},
	{0,207,47},
	{0,205,49},
	{0,203,51},
	{0,201,53},
	{0,199,55},
	{0,197,57},
	{0,195,59},
	{0,193,61},
	{0,191,63},
	{0,189,65},
	{0,187,67},
	{0,185,69},
	{0,183,71},
	{0,181,73},
	{0,179,75},
	{0,177,77},
	{0,175,79},
	{0,173,81},
	{0,171,83},
	{0,169,85},
	{0,167,87},
	{0,165,89},
	{0,163,91},
	{0,161,93},
	{0,159,95},
	{0,157,97},
	{0,155,99},
	{0,153,101},
	{0,151,103},
	{0,149,105},
	{0,147,107},
	{0,145,109},
	{0,143,111},
	{0,141,113},
	{0,139,115},
	{0,137,117},
	{0,135,119},
	{0,133,121},
	{0,131,123},
	{0,129,125},
	{0,127,127},
	{0,125,129},
	{0,123,131},
	{0,121,133},
	{0,119,135},
	{0,117,137},
	{0,115,139},
	{0,113,141},
	{0,111,143},
	{0,109,145},
	{0,107,147},
	{0,105,149},
	{0,103,151},
	{0,101,153},
	{0,99,155},
	{0,97,157},
	{0,95,159},
	{0,93,161},
	{0,91,163},
	{0,89,165},
	{0,87,167},
	{0,85,169},
	{0,83,171},
	{0,81,173},
	{0,79,175},
	{0,77,177},
	{0,75,179},
	{0,73,181},
	{0,71,183},
	{0,69,185},
	{0,67,187},
	{0,65,189},
	{0,63,191},
	{0,61,193},
	{0,59,195},
	{0,57,197},
	{0,55,199},
	{0,53,201},
	{0,51,203},
	{0,49,205},
	{0,47,207},
	{0,45,209},
	{0,43,211},
	{0,41,213},
	{0,39,215},
	{0,37,217},
	{0,35,219},
	{0,33,221},
	{0,31,223},
	{0,29,225},
	{0,27,227},
	{0,25,229},
	{0,23,231},
	{0,21,233},
	{0,19,235},
	{0,17,237},
	{0,15,239},
	{0,13,241},
	{0,11,243},
	{0,9,245},
	{0,7,247},
	{0,5,249},
	{0,3,251},
	{0,1,253}
};


// Function to render particles with coloration based on speed
float current_max=0, last_max = 1000;
float current_min=0, last_min = -1000;
void render_particles_spd(SDL_Renderer* renderer, particle* particles) {
    float denom = (last_max>(-last_min)) ? last_max : -last_min;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            float cif = cbrtf(particles[i].vx/denom);
            int ci = (int)(127.0*cif + 128.0);
            if(ci>255){ ci=255; }
            else if(ci<0){ ci=0; }
            // Colorize based on particle speed
            SDL_SetRenderDrawColor(renderer, color_lut[ci][0], color_lut[ci][1], color_lut[ci][2], 255);
            SDL_Point pts[circle_pixcnt];
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
        }
        if(particles[i].vx > current_max){
            current_max = particles[i].vx;
        }
        if(particles[i].vx < current_min){
            current_min = particles[i].vx;
        }
    }
    last_max = current_max;
    last_min = current_min;
    current_max = 0;
    current_min = 0;
}

// Y speed particle rendering
float current_ymax=0, last_ymax = 1000;
float current_ymin=0, last_ymin = -1000;
void render_particles_yspd(SDL_Renderer* renderer, particle* particles) {
    float denom = (last_ymax>(-last_ymin)) ? last_ymax : -last_ymin;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            float cif = particles[i].vy/denom;
            if(py<HEIGHT/2){ cif = -cif; }
            cif = cbrtf(cif);
            cif = (cif>0) ? 1.0 : -1.0;
            int ci = (int)(127.0*cif + 128.0);
            if(ci>255){ ci=255; }
            else if(ci<0){ ci=0; }
            // Colorize based on particle speed
            SDL_SetRenderDrawColor(renderer, color_lut[ci][0], color_lut[ci][1], color_lut[ci][2], 255);
            SDL_Point pts[circle_pixcnt];
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
        }
        if(particles[i].vx > current_ymax){
            current_ymax = particles[i].vy;
        }
        if(particles[i].vx < current_ymin){
            current_ymin = particles[i].vy;
        }
    }
    //last_ymax = current_ymax;
    //last_ymin = current_ymin;
    current_ymin = 0;
    current_ymax = 0;
}

// Function to render particles with green color
void render_particles_green(SDL_Renderer* renderer, particle* particles) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            if(i==0){
                SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
            }else{
                SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
            }
            SDL_Point pts[circle_pixcnt];
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
        }
    }
}

// Function to render particles with coloration based on pressure
static float maxp=0.0, minp=1.0;
void render_particles_pressure(SDL_Renderer* renderer, particle* particles) {
    float denom = maxp-minp;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            int pind = (int)(px/(RENDER_WIDTH/SECTIONS));
            int ci = (int)(255.0*(pressure[pind]-minp)/denom);
            if(ci<0){ ci=0; }
            ci=255-ci;
            if(i==0){
                SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
            }else{
                SDL_SetRenderDrawColor(renderer, color_lut[ci][0], color_lut[ci][1], color_lut[ci][2], 255);
            }
            SDL_Point pts[circle_pixcnt];
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
        }
    }
}

// Function to render particles with coloration based on density
static float maxd=0.0;
void render_particles_density(SDL_Renderer* renderer, particle* particles) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            int pind = (int)(px/(RENDER_WIDTH/SECTIONS));
            int ci = (int)(255.0*density[pind]/maxd);
            if(ci<0){ ci=0; }
            ci=255-ci;
            if(i==0){
                //ci += 128;
                //ci = (ci>255) ? ci-255 : ci;
                SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
            }else{
                SDL_SetRenderDrawColor(renderer, color_lut[ci][0], color_lut[ci][1], color_lut[ci][2], 255);
            }
            SDL_Point pts[circle_pixcnt];
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
        }
    }
}

// Function to render particles with coloration based on temperature
static float maxt=0.0, mint=0.0;
static float stc = 1.332;

void render_particles_temperature(SDL_Renderer* renderer, particle* particles) {
    float denom = maxt - mint;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            int pind = (int)(px/(RENDER_WIDTH/SECTIONS));
            int ci = (int)(255.0*(temp[pind] - mint)/denom);
            if(ci<0){ ci=0; }
            ci=255-ci;
            if(i==0){
                //ci += 128;
                //ci = (ci>255) ? ci-255 : ci;
                SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
            }else{
                SDL_SetRenderDrawColor(renderer, color_lut[ci][0], color_lut[ci][1], color_lut[ci][2], 255);
            }
            SDL_Point pts[circle_pixcnt];
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
        }
    }
}


void draw_walls(SDL_Renderer* renderer){
    SDL_SetRenderDrawColor(renderer, wall_color[0], wall_color[1], wall_color[2], 255);
    SDL_RenderDrawPoints(renderer, wall_high, wall_cnt_h);
    SDL_RenderDrawPoints(renderer, wall_low, wall_cnt_l);
}

uint64_t getus(){
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    uint64_t tm = (uint64_t)now.tv_sec*1000000 + (uint64_t)now.tv_nsec/1000;
    return tm;
}

static int last_mps = -1;

//main function for rendering and updating simulation
int main(int argc, char* argv[]) {

    setpriority(PRIO_PROCESS, 0, -20);

    //initialize SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("Particle Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, RENDER_WIDTH, HEIGHT, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if(renderer == NULL){ printf("ERR: Failure to create renderer.\n\n"); exit(EXIT_FAILURE); }
    //SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);

    //initialize particles
    particle particles[NUM_PARTICLES];
    init_particles(particles);

    generate_circle(circle_lut, &circle_pixcnt, 0, 0, (int)RADIUS);
    #ifdef DRAW_WALLS
    generate_circle(wall_high, &wall_cnt_h, (int)BIGCIRC_XR, (int)BIGCIRC_YH, (int)WALL_RADIUS);
    generate_circle(wall_low, &wall_cnt_l, (int)BIGCIRC_XR, (int)BIGCIRC_YL, (int)WALL_RADIUS);
    #endif

    // Compute the mean inter-particle distance assuming even distribution
    float ideal_circ_area = (WIDTH*HEIGHT - wall_cnt_h - wall_cnt_l)*0.9069 / (float)NUM_PARTICLES;
    mean_fpl = 2*sqrtf(ideal_circ_area/M_PI) * FPL_MULT;

    printf("Mean free path length is %.3f pixels\n\n", mean_fpl);

    printf("Top big circle center (%f, %f)\n", BIGCIRC_X, BIGCIRC_YH);
    printf("Bottom big circle center (%f, %f)\n\n", BIGCIRC_X, BIGCIRC_YL);
    get_areas();
    printf("\n");

    initial_energy=0.0;

    int pupdate=0;
    int kick=0;
    float slowdown=0;

    int64_t now;
    int64_t tm0 = getus();
    int64_t tm1 = getus();
    int64_t tm_render_particles=0, tm_render_walls=0, tm_present=0, tm_simulate=0, tm_full = 0;
    int tcnt = 0;

    float looptm = 10.0;

    char rc = get_render_color(&bpwr_command);
    int pcnt = 0;

    //main loop for simulation
    while (1) {
        SDL_Event event;
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
            break;
        }

        // Clear screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // Wait for next frame period
        #ifdef FRAMETIME_US
            if(getus() - tm0 > FRAMETIME_US){
                fprintf(stderr, "ERROR: Missed frame timing\n");
            }else{
                while(getus() - tm0 < FRAMETIME_US){
                    usleep(100);
                }
            }
        #endif

        now = getus();
        tm_full += now - tm0;
        tm0 = now;
        tm1 = now;

        // Render particles based on the readut from the start/stop file
        switch (rc){
            case 't':
                render_particles_temperature(renderer, particles);
                break;
            case 'T':
                render_particles_temperature(renderer, particles);
                break;
            case 'p':
                render_particles_pressure(renderer, particles);
                break;
            case 'P':
                render_particles_pressure(renderer, particles);
                break;
            case 'd':
                render_particles_density(renderer, particles);
                break;
            case 'D':
                render_particles_density(renderer, particles);
                break;
            case 's':
                render_particles_spd(renderer, particles);
                break;
            case 'S':
                render_particles_spd(renderer, particles);
                break;
            case 'y':
                render_particles_yspd(renderer, particles);
                break;
            case 'Y':
                render_particles_yspd(renderer, particles);
                break;
            default:
                render_particles_green(renderer, particles);
                break;
        }

        now = getus();
        tm_render_particles += now-tm1;
        tm1 = now;

        // Draw walls
        #ifdef DRAW_WALLS
        draw_walls(renderer);
        #endif
        now = getus();
        tm_render_walls += now-tm1;
        tm1 = now;

        SDL_RenderPresent(renderer);
        now = getus();
        tm_present += now-tm1;
        tm1 = now;

        
        if(kick >= KICK_PERIOD){
            kick=0;
            avg_speed = speed_sum/speed_cnt;
            speed_sum = 0;
            speed_cnt = 0;
            if(init_avg_speed==0){ init_avg_speed = avg_speed; }
            else{
                #ifdef MAGIC_FORCE
                    bpwr = bpwr_command*(init_avg_speed - avg_speed);
                    if(bpwr>0){ bpwr = 0; }
                    for(int i=0; i<NUM_PARTICLES; i++){
                        particles[i].vx += bpwr;
                    }
                #endif
            }
            printf("KICK:\n");
            fflush(stdout);
        }
        kick++;
        
        
        if(pupdate>=PFRAMES){
            pcnt++;
            rc = get_render_color(&bpwr_command);
            pupdate=0;
            maxp = 0.0;
            minp = 1.0;
            maxd = 0.0;
            maxt = 0.0;
            mint = 1.0;
            float mach[SECTIONS];
            float entropy[SECTIONS];
            float min_press = 0;
            int min_press_sec = 0;
            // Adjust each section by scaling factors
            for(int i=0; i<SECTIONS; i++){
                // X speed given per particle
                speeds[i] = x_spds[i]/spd_cnt[i];
                // Avg KE given per particle
                temperature[i] = temperature[i]/spd_cnt[i];
                // Densities given per unit neck area
                spd_cnt[i] = spd_cnt[i]/d_areas[i];
                // Ideal gas law pressure is temp * density
                pressure_igl[i] = temperature[i]*spd_cnt[i];
                // Mach number
                mach[i] = speeds[i]/sqrtf(temperature[i]);
            }
            
            dt = looptm*PPS / (1000.0*(sqrtf(temperature[SECTIONS/2]) - speeds[SECTIONS/2]));

            stc = -1.0*speeds[SECTIONS/2] / sqrtf(temperature[SECTIONS/2]);


            // Get the sums for each vector
            double psum = 0, ssum=0, dsum=0, tsum=0, isum=0, asum=0, msum=0, esum=0;
            for(int i=0; i<SECTIONS; i++){
                psum += bounce_accum[i];
                ssum += speeds[i];
                dsum += spd_cnt[i];
                tsum += temperature[i];
                isum += pressure_igl[i];
                asum += d_areas[i];
                msum += mach[i];
                // Entropy
                entropy[i] = logf(temperature[i]/temperature[SECTIONS-1]) - 0.398*logf(pressure_igl[i]/pressure_igl[SECTIONS-1]);
                esum += entropy[i];
            }
            // Compute relative output vectors and clear
            float spd[SECTIONS];
            for(int i=0; i<SECTIONS; i++){
                pressure[i] = bounce_accum[i]/psum;
                spd[i] = speeds[i]/ssum;
                density[i] = spd_cnt[i]/dsum;
                temp[i] = temperature[i]/tsum;
                pressure_igl[i] = pressure_igl[i]/isum;
                mach[i] = mach[i]*0.84;
                //mach[i] = mach[i]/2.0;
                entropy[i] = entropy[i]/esum;
                bounce_accum[i]=0;
                speeds[i]=0;
                spd_cnt[i]=0;
                temperature[i]=0;
                if(pressure[i]>maxp){ maxp = pressure[i]; }
                if(pressure[i]<minp){ minp = pressure[i]; }
                if(density[i]>maxd){ maxd = density[i]; }
                if(temp[i]>maxt){ maxt = temp[i]; }
                if(temp[i]<mint){ mint = temp[i]; }
            }
            printf("SPEEDS: %d", (int)(spd[0]*1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(spd[i]*1000));
            }
            printf("\n");
            printf("PRESSURES: %d", (int)(pressure[0]*1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(pressure[i]*1000));
            }
            printf("\n");
            printf("DENSITIES: %d", (int)(density[0]*1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(density[i]*1000));
            }
            printf("\n");
            printf("TEMPS: %d", (int)(temp[0]*1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(temp[i]*1000));
            }
            printf("\n");
            printf("IGLP: %d", (int)(pressure_igl[0]*1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(pressure_igl[i]*1000));
            }
            printf("\n");
            printf("AREA: %d", (int)(d_areas[0]/asum*1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(d_areas[i]/asum*1000));
            }
            printf("\n");
            printf("MACH: %d", (int)(mach[0]*-1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(mach[i]*-1000));
            }
            printf("\n");
            printf("ENTROPY: %d", (int)(entropy[0]*1000));
            for(int i=1; i<SECTIONS; i++){
                printf(" %d", (int)(entropy[i]*1000));
            }
            printf("\n");
            looptm = (float)tm_full / 1000.0 / tcnt;
            float simtm = (float)tm_simulate / 1000.0 / tcnt;
            float rentm = ((float)tm_render_particles + (float)tm_render_walls) / 1000.0 / tcnt;
            float prestm = (float)tm_present / 1000.0 / tcnt;
            tcnt=0; tm_full=0; tm_render_particles=0; tm_render_walls=0; tm_simulate=0; tm_present=0;
            printf("SLOWDOWN:%c %.3f, LOOP_TM:%.3fms, ASPD:%.1fppf, STC:%.3fppf, BPWR:%.1fppf, BPWR_COM:%.2fppf, DT:%.4f\n", rc, slowdown, looptm, avg_speed*100, stc, bpwr*100, bpwr_command, dt);
            fflush(stdout);
            // fprintf(stderr, "EBLEED BY SECTION: ");
            // for(int i=0; i<SECTIONS; i++){
            //     fprintf(stderr, "%.2f, ", 1.0-(1.0-temp[i]*SECTIONS)/8.0);
            // }
            // fprintf(stderr, "\n\n");
        }
        pupdate++;

        // Advance the simulation by one timestep
        slowdown = simulate_timestep(particles, mean_fpl);

        now = getus();
        tm_simulate += now-tm1;
        tm1 = now;
        tcnt++;
        //usleep(3000);
    }

    //clean up SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}