#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL2/SDL.h>
#include <sys/resource.h>
#include <time.h>

const char control_fname[] = "/home/jtmcg/psim/render_color.txt";

// Window geometry
#define WIDTH 3680.0
#define RENDER_WIDTH 3680.0
#define HEIGHT 1050.0
#define WALL_RADIUS 3500.0
#define LITTLE_RADIUS (WALL_RADIUS/10.0)
#define NECK_HEIGHT 400.0

// Simulation parameters
#define NUM_PARTICLES 10000 //number of particles in the simulation
#define RESERVE_PARTICLES 5000
#define RADIUS 3.5 //radius of each particle
#define DT 0.005 //timestep for simulation
#define REPOS 1.00001 // Repositioning constant. Puts particles slightly further away on bounce to prevent them getting stuck

// Initial temperature, scales the randomized particle starting speeds
#define INITIAL_TEMP 650.0
#define PRESSURE_RATIO 0.5
// Blower power, expressed as a fraction of initial temperature

// Surface roughness, scales the random shifts applied to wall bounce angle
#define ROUGHNESS_ANG 0.0*M_PI/180.0
// Number of distinct sections for pressure and speed computations
#define SECTIONS 24
#define PFRAMES 400
#define KICK_PERIOD 80

#define FORCE 0.0
#define FPL_MULT 0.8

//#define INIT_HALF

// Define this to draw walls. Takes a while on startup to compute the draw circle
//#define DRAW_WALLS
const int wall_color[3] = {255, 0, 0};

// Define the frame time target in usecs
//#define FRAMETIME_US 20833

#define ROOT2O2 0.7071067811865476



//structure for a single particle
typedef struct {
    float x; //x position of particle
    float y; //y position of particle
    float vx; //x velocity of particle
    float vy; //y velocity of particle
} particle;

static particle particles[NUM_PARTICLES];
static int reserve_inds[NUM_PARTICLES];
static int reserve_count = RESERVE_PARTICLES;

// Returns a random float between 0 and 1
inline float randf(){ return ((float)rand() / RAND_MAX); }

// Function to read single character from control file
char get_render_color(){
    FILE* fp = fopen(control_fname, "r");
    if( fp == NULL ){
        return 0;
    }
    int char0 = fgetc(fp);
    if( char0==EOF || char0=='\n' ){
        fclose(fp);
        return 0;
    }
    fclose(fp);
    return (char)char0;
}

// Function to randomize the input bounce angle according to the roughness angle
inline void rand_bounce(float *nx, float *ny){
    float rang = (randf()-0.5)*2*ROUGHNESS_ANG;
    float c = cosf(rang);
    float s = sinf(rang);
    float ox = *nx;
    float oy = *ny;
    *nx = ox*c - oy*s;
    *ny = ox*s + oy*c;
}

// Householder reflection function for particle velocities
inline void householder(particle *p, float nx, float ny){
    float s = 2*(nx*p->vx + ny*p->vy);
    p->vx -= nx*s;
    p->vy -= ny*s;
}

// Geometry definitions for walls
#define BIGCIRC_X (WIDTH/2)
#define UNREND ((WIDTH-RENDER_WIDTH)/2.0)
#define BIGCIRC_YH (HEIGHT/2+NECK_HEIGHT/2+WALL_RADIUS)
#define BIGCIRC_YL (HEIGHT/2-NECK_HEIGHT/2-WALL_RADIUS)
#define PUSH_RAD ((WALL_RADIUS+RADIUS)*REPOS)

// Measurement arrays
static float bounce_accum[SECTIONS];
static float speeds[SECTIONS];
static float last_speeds[SECTIONS];
static float pressure[SECTIONS];
static float density[SECTIONS];
static float temp[SECTIONS];
static float temperature[SECTIONS];
static float spd_cnt[SECTIONS];
static float p_areas[SECTIONS];
static float d_areas[SECTIONS];

static float mean_fpl;
static float circle_intersect;

// Function to fill the "areas" for relative pressure and density computations
#define BIGCIRC_XR (BIGCIRC_X-UNREND)
void get_areas(){
    float lb, rb, cx, cy, ar, nar;
    circle_intersect = sqrt(WALL_RADIUS*WALL_RADIUS - BIGCIRC_YL*BIGCIRC_YL);
    float theta0 = asinf(circle_intersect/WALL_RADIUS);
    float theta1 = asinf(-BIGCIRC_YL/WALL_RADIUS);
    float lc_bound = BIGCIRC_XR - circle_intersect;
    float rc_bound = BIGCIRC_XR + circle_intersect;
    // Iterate through the sections
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
            // Length along wall plus arc length of included circle section
            ar = (lc_bound-lb) + WALL_RADIUS*(theta0 - theta);
            // This is the trapezoidal approximation for average area
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
        // Save the pressure areas
        p_areas[i] = 2*ar;
        // Save the average channel heights used for density computation
        d_areas[i] = nar;
        printf("SECTION %d: WALL AREA %.1f, NECK AREA %.1f\n", i, 2*ar, nar);
    }
}

static float feed_rate = 0;
static float bleed_rate = 0;
static double upstream_psum = 0;
static double downstream_psum = 0;
static int chance_cnt = 0;

// Function to execute bounces for particles striking channel walls
uint8_t venturi_bounce(particle *p, int i){
    // Gather information for statistics
    float px = p->x;
    int pind = (int)((px - UNREND)/(RENDER_WIDTH/SECTIONS));
    uint8_t pvalid = (px < RENDER_WIDTH+UNREND) && (px >= UNREND);
    float vx = p->vx, vy = p->vy;
    if(pvalid){
        speeds[pind] += vx;
        vx -= last_speeds[pind];
        temperature[pind] += (vx*vx + vy*vy)/8.0;
        spd_cnt[pind] += 1.0;
    }

    float s, nx, ny;
    uint8_t bounds_exceed = 0;

    // Left edge bounce/dissapear behavior
    if(px <= RADIUS){
        bounds_exceed = 1;
        chance_cnt++;
        // Randomly determine if this particle will be released
        if(randf() < bleed_rate){
            // Was released
            // Add index to reserve list and nullify the particle params
            if(reserve_count<NUM_PARTICLES){
                reserve_inds[reserve_count] = i;
                reserve_count++;
            }else{
                fprintf(stderr, "Last particle was bled. Something is very wrong...\n\n");
                exit(1);
            }
            p->x = -1;
            p->y = -1;
            p->vx = 0;
            p->vy = 0;
        // Otherwise normal bounce
        }else{
            downstream_psum = (p->vx < 0) ? (downstream_psum - p->vx) : (downstream_psum + p->vx);
            nx=-1;
            ny=0;
            rand_bounce(&nx, &ny);
            householder(p, nx, ny);
            p->x = RADIUS*REPOS;
        }
    }
    // Right edge bounce behavior
    else if(px >= (WIDTH-RADIUS)){
        bounds_exceed = 1;
        upstream_psum = (p->vx < 0) ? (upstream_psum - p->vx) : (upstream_psum + p->vx);
        nx=1;
        ny=0;
        rand_bounce(&nx, &ny);
        householder(p, nx, ny);
        p->x = WIDTH-RADIUS*REPOS;
    }

    float py = p->y;
    // Check Y rectangle bounds
    if(py <= RADIUS){
        bounds_exceed = 1;
        float ay = p->vy;
        float ax = p->vx;
        // Increment the correct bounce accumulator for pressure graph
        if(pvalid){
            bounce_accum[pind] -= ay;
        }
        // Increment the upstream/downstream pressure accumulators if x in range
        // float px = p->x;
        // if(px<(BIGCIRC_X-circle_intersect)){
        //     downstream_psum = (ay < 0) ? (downstream_psum - ay) : (downstream_psum + ay);
        // }else if(px>(BIGCIRC_X+circle_intersect)){
        //     upstream_psum = (ay < 0) ? (upstream_psum - ay) : (upstream_psum + ay);
        // }
        nx=0;
        ny=-1;
        rand_bounce(&nx, &ny);
        householder(p, nx, ny);
        p->y = RADIUS*REPOS;
    }
    else if(py >= (HEIGHT-RADIUS)){
        bounds_exceed = 1;
        float ay = p->vy;
        float ax = p->vx;
        if(pvalid){
            bounce_accum[pind] += ay;
        }
        // Increment the upstream/downstream pressure accumulators if x in range
        // float px = p->x;
        // if(px<(BIGCIRC_X-circle_intersect)){
        //     downstream_psum = (ay < 0) ? (downstream_psum - ay) : (downstream_psum + ay);
        // }else if(px>(BIGCIRC_X+circle_intersect)){
        //     upstream_psum = (ay < 0) ? (upstream_psum - ay) : (upstream_psum + ay);
        // }
        nx=0;
        ny=1;
        rand_bounce(&nx, &ny);
        householder(p, nx, ny);
        p->y = HEIGHT-RADIUS*REPOS;
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
        if(mom<0){ mom=0-mom; }
        bounce_accum[pind] += mom;
        p->x = BIGCIRC_X + PUSH_RAD*nx;
        p->y = BIGCIRC_YH + PUSH_RAD*ny;
        rand_bounce(&nx, &ny);
        householder(p, nx, ny);
        //householder(p, nx, ny, 1.0);
        return 1;
    }

    dy = p->y - BIGCIRC_YL;
    dist = sqrt(dx2 + dy*dy);
    if( dist<(WALL_RADIUS+RADIUS) ){
        // Perform householder reflection of velocity
        nx = dx/dist;
        ny = dy/dist;
        mom = nx*p->vx + ny*p->vy;
        if(mom<0){ mom=0-mom; }
        bounce_accum[pind] += mom;
        p->x = BIGCIRC_X + PUSH_RAD*nx;
        p->y = BIGCIRC_YL + PUSH_RAD*ny;
        rand_bounce(&nx, &ny);
        householder(p, nx, ny);
        //householder(p, nx, ny, 1.0);
        return 1;
    }

    return bounds_exceed;
}

// Function to check if particles overlap with themselves or bounds
inline uint8_t check_overlap(float x, float y, int occ){
    particle p;
    p.x=x; p.y=y; p.vx=0; p.vy=0;
    #ifdef INIT_HALF
    if(x<WIDTH*0.65 || x>WIDTH*0.9){ return 1; }
    #endif
    if(venturi_bounce(&p, 0)){ return 1; }
    for(int i=0; i<occ; i++){
        float dx = particles[i].x - x;
        float dy = particles[i].y - y;
        float dist = sqrt(dx*dx + dy*dy);
        if(dist<=2*RADIUS){ return 1; }
    }
    return 0;
}

// Function to initialize particles with random positions and velocities
void init_particles() {
    int i=0;
    while(i<(RESERVE_PARTICLES)){
        reserve_inds[i]= i;
        particles[i].x = -1;
        particles[i].y = -1;
        particles[i].vx = 0;
        particles[i].vy = 0;
        i++;
    }
    i = RESERVE_PARTICLES;
    reserve_count = RESERVE_PARTICLES;
    while(i<NUM_PARTICLES){
        reserve_inds[i] = -1;
        float cx = randf()*WIDTH;
        float cy = randf()*HEIGHT;
        if(!check_overlap(cx, cy, i)){
            particles[i].x = cx;
            particles[i].y = cy;
            particles[i].vy = (randf()*2 - 1)*INITIAL_TEMP;
            particles[i].vx = (randf()*2 - 1)*INITIAL_TEMP;
            i++;
        }
    }
}


// Function to simulate a timestep for each particle
float simulate_timestep() {

    // Scan through all particle pair combos
    for (int i = 0; i < NUM_PARTICLES - 1; i++) {
        float ix = particles[i].x;
        float iy = particles[i].y;
        if(iy<0 && ix<0){ continue; }
        float ivx = particles[i].vx;
        float ivy = particles[i].vy;
        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            float jx = particles[j].x;
            float jy = particles[j].y;
            if(iy<0 && ix<0){ continue; }
            float jvx = particles[j].vx;
            float jvy = particles[j].vy;
            // Compute distance with distance formula
            float dx = ix - jx;
            float dy = iy - jy;
            float dist2 = dx*dx + dy*dy;
            float dist = sqrtf(dist2);
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
                    float force = dist*FORCE;
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

    // Update positions of particles
    for (int i = 0; i < NUM_PARTICLES; i++) {
        if(particles[i].x<0 && particles[i].y<0){ continue; }
        float vx = particles[i].vx;
        float vy = particles[i].vy;
        particles[i].x += vx * DT;
        particles[i].y += vy * DT;
        // Execute wall bounces
        venturi_bounce(particles+i, 1.0);
    }

    // Return the current amount of slowdown to maintain energy
    return 0;
}

// Define a table to hold pixel offsets for rendering a filled in circle
static SDL_Point circle_lut[(int)((RADIUS+1)*(RADIUS+1)*4)];
static int circle_pixcnt=0;

static SDL_Point wall_high[(int)(WALL_RADIUS*WALL_RADIUS*M_PI/2.0)];
static int wall_cnt_h=0;
static SDL_Point wall_low[(int)(WALL_RADIUS*WALL_RADIUS*M_PI/2.0)];
static int wall_cnt_l=0;

// Function to generate tables of pixels for filled in circles
void generate_circle(SDL_Point *lut, int *pixcnt, int x_off, int y_off, int rad){
    float ang=0.0;
    float xf, yf;
    int x, y;
    uint8_t use_bounds = (x_off!=0 || y_off!=0);
    int radc = (rad*11)/10;
    
    float ang_inc = 0.8 / (float)rad;
    printf("\nUSE BOUNDS IS %d, ANGLE INCREMENT IS %.6f RAD\n\n", (int)use_bounds, ang_inc);
    while(ang<(2*M_PI)){
        xf = cosf(ang);
        yf = sinf(ang);
        ang += ang_inc;
        //printf("Angle is %.3f of %.3f, pixcnt %d\n", ang, 2*M_PI, *pixcnt);
        for(float rf=0; rf<=rad; rf+=0.65){
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
                if(!use_bounds || (x>=0 && x<RENDER_WIDTH && y>=0 && y<HEIGHT)){
                    lut[(*pixcnt)].x = x;
                    lut[(*pixcnt)].y = y;
                    (*pixcnt) ++;
                }
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
void render_particles_spd(SDL_Renderer* renderer) {
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
    //last_max = current_max;
    //last_min = current_min;
    current_max = 0;
    current_min = 0;
}

// Y speed particle rendering
float current_ymax=0, last_ymax = 1000;
float current_ymin=0, last_ymin = -1000;
void render_particles_yspd(SDL_Renderer* renderer) {
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
void render_particles_green(SDL_Renderer* renderer) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
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
void render_particles_pressure(SDL_Renderer* renderer) {
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

// Function to render particles with coloration based on density
static float maxd=0.0;
void render_particles_density(SDL_Renderer* renderer) {
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
void render_particles_temperature(SDL_Renderer* renderer) {
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

void feed_particles(int cnt, float fac){
    uint8_t warn = 0;
    for(int i=0; i<cnt; i++){
        reserve_count--;
        if(reserve_count<0){
            reserve_count = 0;
            warn = 1;
        }else{
            int ind = reserve_inds[reserve_count];
            particles[ind].x = WIDTH-(7*RADIUS);
            particles[ind].y = 2*RADIUS + (HEIGHT-4*RADIUS)*randf();
            particles[ind].vy = (randf()*2 - 1)*INITIAL_TEMP*fac;
            particles[ind].vx = (randf()*2 - 1)*INITIAL_TEMP*fac;
        }
    }
    if(warn){
        fprintf(stderr, "WARNING: RAN OUT OF RESERVE PARTICLES. BOOST RESERVE COUNT\n\n");
    }
}

//main function for rendering and updating simulation
int main(int argc, char* argv[]) {

    setpriority(PRIO_PROCESS, 0, -20);

    //initialize SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("Particle Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, RENDER_WIDTH, HEIGHT, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    //initialize particles
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

    int pupdate=0;
    int kick=0;
    float slowdown=0;

    int64_t now;
    int64_t tm0 = getus();
    int64_t tm1 = getus();
    int64_t tm_render_particles=0, tm_render_walls=0, tm_present=0, tm_simulate=0, tm_full = 0;
    int tcnt = 0;

    char rc = get_render_color();

    double initial_pressure = 0;
    float initial_temp = 0, ftemp=0;
    float feed_fac = 0;
    float pressure_ratio = 1.0;
    int feed_cnt = 0;
    float upstream_pressure, downstream_pressure;

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
                render_particles_temperature(renderer);
                break;
            case 'T':
                render_particles_temperature(renderer);
                break;
            case 'p':
                render_particles_pressure(renderer);
                break;
            case 'P':
                render_particles_pressure(renderer);
                break;
            case 'd':
                render_particles_density(renderer);
                break;
            case 'D':
                render_particles_density(renderer);
                break;
            case 's':
                render_particles_spd(renderer);
                break;
            case 'S':
                render_particles_spd(renderer);
                break;
            case 'y':
                render_particles_yspd(renderer);
                break;
            case 'Y':
                render_particles_yspd(renderer);
                break;
            default:
                render_particles_green(renderer);
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
            if(initial_pressure==0){
                initial_pressure = (upstream_psum+downstream_psum)/2;
                initial_temp = temp[SECTIONS-1];
                feed_fac = chance_cnt*4;
            }else{
                // Compute relative upstream and downstream pressures for logging
                upstream_pressure = upstream_psum/initial_pressure;
                downstream_pressure = downstream_psum/initial_pressure;
                //pressure_ratio = downstream_psum/upstream_psum;
                pressure_ratio = downstream_psum/upstream_psum;
                // Compute the particle bleed rate (probability of particle exit)
                bleed_rate = (downstream_pressure-PRESSURE_RATIO);
                // Compute the number of particles to inject
                feed_cnt = (int)(feed_fac*(1.0-upstream_pressure));
                if(feed_cnt>RESERVE_PARTICLES-2){ feed_cnt = RESERVE_PARTICLES-2; }
                if(initial_temp>0){
                    feed_particles(feed_cnt, initial_temp/ftemp*2.0);
                }
            }
            
            // Reset accumulators for next go around
            upstream_psum = 0;
            downstream_psum = 0;
            chance_cnt = 0;
            printf("KICK:\n");
        }
        kick++;
        
        if(pupdate>=PFRAMES){
            rc = get_render_color();
            pupdate=0;
            maxp = 0.0;
            minp = 1.0;
            maxd = 0.0;
            maxt = 0.0;
            mint = 1.0;
            // Adjust each section by scaling factors
            for(int i=0; i<SECTIONS; i++){
                // X speed given per particle
                speeds[i] = speeds[i]/spd_cnt[i];
                last_speeds[i] = speeds[i];
                // Temp given per particle
                temperature[i] = temperature[i]/spd_cnt[i];
                // Bounce accumulation given per unit wall area
                bounce_accum[i] = bounce_accum[i]/p_areas[i];
                // Densities given per unit neck area
                spd_cnt[i] = spd_cnt[i]/d_areas[i];
            }

            if(initial_temp==0){
                initial_temp = temperature[SECTIONS-1];
            }
            ftemp = temperature[SECTIONS-1];
            
            // Get the sums for each vector
            double psum = 0, ssum=0, dsum=0, tsum=0;
            for(int i=0; i<SECTIONS; i++){
                psum += bounce_accum[i];
                ssum += speeds[i];
                dsum += spd_cnt[i];
                tsum += temperature[i];
            }
            // Compute relative output vectors and clear
            float spd[SECTIONS];
            for(int i=0; i<SECTIONS; i++){
                pressure[i] = bounce_accum[i]/psum;
                spd[i] = speeds[i]/ssum;
                density[i] = spd_cnt[i]/dsum;
                temp[i] = temperature[i]/tsum;
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
            float looptm = (float)tm_full / 1000.0 / tcnt;
            float simtm = (float)tm_simulate / 1000.0 / tcnt;
            float rentm = ((float)tm_render_particles + (float)tm_render_walls) / 1000.0 / tcnt;
            float prestm = (float)tm_present / 1000.0 / tcnt;
            tcnt=0; tm_full=0; tm_render_particles=0; tm_render_walls=0; tm_simulate=0; tm_present=0;
            printf("SLOWDOWN:%c %.3f, LOOP_TM:%.3fms, SIM_TM:%.3fms, UP_P:%3f, DWN_P:%.3f, PR:%.3f, BLEED:%.3f, FEED_CNT:%d\n", rc, initial_temp/ftemp, looptm, simtm, upstream_pressure, downstream_pressure, pressure_ratio, bleed_rate, (int)feed_cnt);
            fflush(stdout);
            // fprintf(stderr, "EBLEED BY SECTION: ");
            // for(int i=0; i<SECTIONS; i++){
            //     fprintf(stderr, "%.2f, ", 1.0-(1.0-temp[i]*SECTIONS)/8.0);
            // }
            // fprintf(stderr, "\n\n");
        }
        pupdate++;

        // Advance the simulation by one timestep
        slowdown = simulate_timestep();

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