#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <sys/resource.h>
#include <time.h>

const char control_fname[] = "/home/jtmcg/projects/psim/ctrl.txt";

// Window geometry
#define WIDTH 4500.0
#define RENDER_WIDTH 3680.0
#define HEIGHT 900.0
#define UNREND ((WIDTH-RENDER_WIDTH)/2)
#define TXT_HEIGHT 45
// Channel geometry
#define THROAT_WIDTH (HEIGHT*0.45)
//#define THROAT_WIDTH  HEIGHT
#define CONVERGE_START 0.9
#define THROAT_POS 0.25
#define DIVERGE_START 0.02
#define DIVERGE_END 0.01

// Simulation parameters
#define VSYNC_PARAM -1     // Change VSYNC param.  1 = unrestricted,  -1 = vsync,  0 = ???
#define NUM_PARTICLES 2800 // number of particles in the simulation
#define EXIT_PX 0.85
#define RADIUS 6.7   // radius of each particle
#define INIT_DT 0.18 // initial timestep for simulation
#define REPOS 1.0001 // Repositioning constant. Puts particles slightly further away on bounce to prevent them getting stuck

// Blower power, expressed as a fraction of initial temperature
#define SMOOTHY 0.1

// Surface roughness, scales the random shifts applied to wall bounce angle
//     Effectively generates entropy
#define ROUGHNESS_ANG 0.0

// Number of distinct sections for pressure and speed computations
#define SECTIONS 21
#define PFRAMES 230
#define KICK_PERIOD 20

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

// Function to read information from control file
char get_file_info(float *sim_spd, float *expx){
    FILE* fp = fopen(control_fname, "r");
    if( fp == NULL )  return 0;
    // Pull first line
    char rs[40];
    if(fgets(rs, 40, fp)==NULL)  return 0;
    // Capture the last alpha char value in the string and replace all with white space
    char charr = '\0';
    for(int i=0; i<40; i++){
        char c = rs[i];
        if(c==0)  break;
        else if(c>=65 && c<=90){  charr = c;  rs[i]=' ';  }
        else if(c>=97 && c<=122){  charr = c;  rs[i]=' ';  }
    }
    // Assign value if valid float found in first line
    float ret = atof(rs);
    if(ret!=0.0) *sim_spd = ret;
    // Pull the second line and assign if valid
    if(fgets(rs, 40, fp)==NULL)  return charr;
    ret = atof(rs);
    *expx = ret;

    fclose(fp);
    printf("C:%c SPD:%.2f exPX:%.2f\n", charr, *sim_spd, *expx);
    return charr;
}

// Function to determine which section the particle is in. Negative value is out of bounds
inline int get_section(particle *p){
    float s = floorf((p->x - UNREND)/(RENDER_WIDTH/SECTIONS));
    int si = (int)s;
    if(si>(SECTIONS-1))  return -1;
    return si;
}

static float reflect_rms = 0.0;
void inline lambert_reflect(particle *p, float surf_nx, float surf_ny, float bleed){
    // Generate random lambertian normal vector
    float rf = randf();
    float r1r = rf*(1.0-rf);
    float nx = 2.0*sqrtf(r1r)*(2.0*rf-1.0);
    float ny = 4.0*r1r;
    // Rotate lambertian random normal to surface alignement
    float nx_rot = surf_ny*nx - surf_nx*ny;
    float ny_rot = surf_nx*nx + surf_ny*ny;
    // Keep a running average RMS impact velocity
    float vx = p->vx, vy=p->vy;
    reflect_rms = 0.9*reflect_rms + 0.1*sqrtf(vx*vx + vy*vy);
    // Diffuse reflect the particle
    p->vx = bleed*nx_rot*reflect_rms;
    p->vy = bleed*ny_rot*reflect_rms;
}

// Householder reflection function for particle velocities
static double momentum_sums[SECTIONS] = {0};
inline void householder(particle *p, float nx, float ny){
    // Compute velocity dot product with normal
    float dot = (nx*p->vx + ny*p->vy);
    // Compute section summations for wall impact momentum
    int sect = get_section(p);
    if(sect >= 0)  momentum_sums[sect] += fabsf(dot);
    // if(randf()>0.7){
    //     lambert_reflect(p, nx, ny, bleed);
    //     return;
    // }
    // Set bounce velocity
    float s = 2*dot;
    p->vx -= nx*s;
    p->vy -= ny*s;
}

inline uint8_t top_bot_bounce(particle *p, const float wall_offset){
    float py = p->y;
    // Check Y rectangle bounds
    if(py <= (wall_offset + RADIUS)){
        householder(p, 0, -1);
        p->y = RADIUS*REPOS + wall_offset;
        return 1;
    }
    if(py >= (HEIGHT-RADIUS-wall_offset)){
        householder(p, 0, 1);
        p->y = HEIGHT - RADIUS*REPOS - wall_offset;
        return 1;
    }
    return 0;
}

inline uint8_t constriction(particle *p, const uint8_t open_close, const float throat_width, const float length, const float throat_x){
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
            householder(p, -nx, ny);
            p->y = yc+(REPOS-1.0);
            return 1;
        }
        // Repeat for lower boundary
        const float b_dwn = (-m)*(length-throat_x) + HEIGHT - RADIUS;
        yc = (-m)*(p->x) + b_dwn;
        if(p->y > yc){
            householder(p, -nx, -ny);
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
            householder(p, nx, ny);
            p->y = yc+(REPOS-1.0);
            return 1;
        }
        // Repeat for lower boundary
        const float b_dwn = (-m)*throat_x + b + throat_width - RADIUS;
        yc = m*(p->x) + b_dwn;
        if(p->y > yc){
            householder(p, nx, -ny);
            p->y = yc-(REPOS-1.0);
            return 1;
        }
    }
    return 0;
}

static double speed_sum = 0;
static int64_t speed_cnt = 0;
static float vx_sums[SECTIONS] = {0};
static float vx_means[SECTIONS] = {0};
static float e_sums[SECTIONS] = {0};
static uint32_t pcnts[SECTIONS] = {0};
static float section_areas[SECTIONS] = {0};
static int mass_flow_sum = 0;
static uint32_t mass_flow_cnt = 0;

static float blow_pwr = 0;
uint8_t venturi_bounce(particle *p){
    float x = p->x;
    if(( x < 0 ) || ( x > WIDTH )){
        p->x = -1.0;
        p->y = -1.0;
        p->vx = 0;
        p->vy = 0;
        return 1;
    }
    int s = get_section(p);
    float vx = p->vx, vy = p->vy;
    speed_sum += vx;
    speed_cnt ++;
    // Valid section
    if(s >= 0){
        vx_sums[s] += vx;
        // Subtract bulk particle speed for temp computation
        float vxd = vx - vx_means[s];
        // Increment particle enrgy summation for section
        e_sums[s] += (vxd*vxd + vy*vy);
        pcnts[s] ++;
    }
    // Compute the X offset
    x = p->x - UNREND;
    // Wide wall bounce
    if( x < RENDER_WIDTH*DIVERGE_END ){
        return top_bot_bounce(p, 0.0);
    // Divergent section
    }else if( x < RENDER_WIDTH*DIVERGE_START ){
        return constriction(p, 1, THROAT_WIDTH, RENDER_WIDTH*(DIVERGE_START-DIVERGE_END), RENDER_WIDTH*DIVERGE_START+UNREND);
    // Throat
    }else if( x < RENDER_WIDTH*THROAT_POS ){
        return top_bot_bounce(p, (HEIGHT-THROAT_WIDTH)*0.5);
    // Convergent section
    }else if( x < RENDER_WIDTH*CONVERGE_START ){
        return constriction(p, 0, THROAT_WIDTH, RENDER_WIDTH*(CONVERGE_START-THROAT_POS), RENDER_WIDTH*THROAT_POS+UNREND);
    // Other wide wall
    }else{
        return top_bot_bounce(p, 0.0);
    }
}

// Function to check if particles overlap with themselves or bounds
inline uint8_t check_overlap(particle* particles, float x, float y, int occ){
    particle p;
    p.x=x; p.y=y; p.vx=0; p.vy=0;
    // Check if particle outside of bounds
    if(venturi_bounce(&p)){ return 1; }
    // Check if particle overlaps another
    for(int i=0; i<occ; i++){
        float dx = particles[i].x - x;
        float dy = particles[i].y - y;
        float dist = sqrt(dx*dx + dy*dy);
        if(dist<=2*RADIUS){ return 1; }
    }
    return 0;
}

static float dt = INIT_DT;
// Function to initialize particles with random positions and velocities
void init_particles(particle* particles, float *areas) {
    int i=0;
    // Continue to generate new particles at random locations until desired quantity is reached
    while(i<NUM_PARTICLES){
        float cx = randf()*WIDTH;
        float cy = randf()*HEIGHT;
        if(!check_overlap(particles, cx, cy, i)){
            particles[i].x = cx;
            particles[i].y = cy;
            particles[i].vy = (randf()*2 - 1);
            particles[i].vx = (randf()*2 - 1);
            i++;
        }
    }
    while(i<NUM_PARTICLES){
        particles[i].x = -1.0;
        particles[i].y = -1.0;
        particles[i].vy = 0;
        particles[i].vx = 0;
        i++;
    }
}

inline int null_particle(particle *p){
    float x = p->x;
    if(isnan(x) || x<0){
        p->x = -1.0;
        p->y = -1.0;
        p->vx = 0;
        p->vy = 0;
        return 1;
    }
    return 0;
}

// Declare and zero the initial energy
static float last_dist[NUM_PARTICLES*NUM_PARTICLES/2];
inline void interparticle_bounce(particle *particles){
    // Scan through all particle pair combos
    for (int i = 0; i < (NUM_PARTICLES - 1); i++) {
        float ix = particles[i].x;
        if(ix<0)  continue;
        float iy = particles[i].y;
        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            float jx = particles[j].x;
            if(jx<0)  continue;
            float jy = particles[j].y;
            // Compute distance with distance formula
            float dx = ix - jx;
            float dy = iy - jy;
            float dist2 = dx*dx + dy*dy;
            // Check for collision
            if (dist2 <= (4*RADIUS*RADIUS)) {
                // Get velocities
                float jvx = particles[j].vx;
                float jvy = particles[j].vy;
                float ivx = particles[i].vx;
                float ivy = particles[i].vy;
                // Compute normal vector components
                float dist = sqrt(dist2);
                float nx = dx/dist;
                float ny = dy/dist;
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
            }
        }
    }
}

// Function to return the initial relative pressure
void inline get_num_temp(particle *particles, int *cnt, double *temp){
    double t = 0;
    int c = 0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        float x = particles[i].x;
        if((x >= 0) && ((x<UNREND) || (x>(UNREND+RENDER_WIDTH)))){
            float spd2 = particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy;
            t += spd2;
            c++;
        }
    }
    // Return the average temp/cnt for each half
    *temp = t / c;
    *cnt = c / 2;
}

inline int inject_particles(particle *particles, float injection_energy, int cnt, const uint8_t entry_true){
    // Return null if all requested particles injected
    if(cnt<=0)  return 0;
    // Scan for empty particles
    for(int i=0; i<NUM_PARTICLES; i++){
        float x = particles[i].x;
        float y = particles[i].y;
        if(x<0 && y<0){
            // Inject particle to entry or exit
            particles[i].x = entry_true ? (WIDTH-RADIUS) : RADIUS;
            // Randomize Y position
            particles[i].y = RADIUS + (HEIGHT-2*RADIUS)*randf();
            // Provide random initial Y velocity
            particles[i].vy = injection_energy*(randf()*2.0 - 1.0);
            // Provide initial X velocity directed to the middle
            particles[i].vx = injection_energy*randf()*(entry_true ? -1.0 : 1.0);
            // Tally particle injected
            cnt --;
        }
        // Return null if all requested particles injected
        if(cnt<=0)  return 0;
    }
    // Otherwise return for not enough particles available
    return -1;
}

static double init_temp=0.0;
static int init_cnt = 0;
static int out_cnt_targ = 0;
static float entry_inject_temp = 1.0;
static float exit_inject_temp = 1.0;
static int entry_cnt = 0;
static int exit_cnt = 0;
static float exit_px = EXIT_PX;
// Function to simulate a timestep for each particle
float simulate_timestep(particle* particles, double elapsed_time_ms){
    // Get the initial temperature and density for inlet/outlet regions
    if(init_cnt==0){
        get_num_temp(particles, &init_cnt, &init_temp);
    }
    out_cnt_targ = init_cnt*exit_px;
    // Check all particle pairs for collision events
    interparticle_bounce(particles);
    entry_cnt = 0;
    exit_cnt = 0;
    double entry_energy = 0.0;
    double exit_energy = 0.0; 
    // Update positions of particles
    for (int i = 0; i < NUM_PARTICLES; i++) {
        // Execute wall bounces and LR edge particle removal
        venturi_bounce(particles+i);
        // Filter for present particles
        float x = particles[i].x;
        if(isnan(x)){
            particles[i].x = -1.0;
            particles[i].y = -1.0;
            particles[i].vx = 0;
            particles[i].vy = 0;
        }else if(x >= 0 ){
            // Propagate particle position from speed
            float vx = particles[i].vx;
            float vy = particles[i].vy;
            particles[i].x += (vx * dt * elapsed_time_ms);
            particles[i].y += (vy * dt * elapsed_time_ms);
            // Compute various energies
            if(x<UNREND){
                float spd2 = vx*vx + vy*vy;
                exit_energy += spd2;
                exit_cnt++;
            }else if(x>(RENDER_WIDTH+UNREND)){
                float spd2 = vx*vx + vy*vy;
                entry_energy += spd2;
                entry_cnt++;
            }
        }
    }

    // Adjust injection temperatures to maintain initial state
    if(entry_energy > 0)  entry_inject_temp *= (init_temp*entry_cnt/entry_energy);
    if(exit_energy > 0)  exit_inject_temp *= (init_temp*exit_cnt/exit_energy);
    if(entry_inject_temp > 5)  entry_inject_temp = 5;
    if(exit_inject_temp > 5)  exit_inject_temp = 5;

    // Inject particles to exit
    inject_particles(particles, exit_inject_temp, (out_cnt_targ - exit_cnt), 0);
    // Inject particles to entry
    inject_particles(particles, entry_inject_temp, (init_cnt - entry_cnt), 1);
}

// Define a table to hold pixel offsets for rendering a filled in circle
static SDL_Point circle_lut[(int)((RADIUS+2)*(RADIUS+2)*4)];
static int circle_pixcnt=0;

// Function to generate tables of pixels for filled in circles
void generate_circle(SDL_Point *lut, int *pixcnt, int rad){
    int pcnt = 0;
    float ang=0.0;
    float xf, yf;
    int x, y;
    float ang_inc = 0.3 / (float)rad;
    int radc = rad*2;
    while(ang<(2*M_PI)){
        xf = cosf(ang);
        yf = sinf(ang);
        ang += ang_inc;
        for(float rf=0; rf<=rad; rf+=0.35){
            x = roundf(xf*rf);
            y = roundf(yf*rf);
            // Set flag if pixel is already accounted for
            uint8_t present=0;
            for(int i=0; i<pcnt; i++){
                if((lut[i].x == x) && (lut[i].y == y)){  present=1;  break;  }
            }
            // Place pixel into LUT if not present
            if(!present){
                lut[pcnt].x = x;
                lut[pcnt].y = y;
                pcnt++;
            }
        }
    }
    *pixcnt = pcnt;
}


static uint32_t wl_actual = 0;
static SDL_Point wall_lut[(int)(RENDER_WIDTH*HEIGHT)];
// Function to build a table of points to render for the walls
void build_wall(float *areas){
    // Zero section counts
    uint32_t section_null_cnts[SECTIONS];
    memset(section_null_cnts, 0, sizeof(section_null_cnts));
    // Zero the wall LUT
    memset(wall_lut, 0, sizeof(wall_lut));
    uint32_t walli = 0;
    particle p;
    p.vx=0; p.vy=0;
    // Iterate through all locations in the renderspace
    for(int x=0; x<RENDER_WIDTH; x++){
        p.x=x+UNREND;
        for(int y=0; y<HEIGHT; y++){
            p.y=y;
            // Draw point if particle outside of bounds
            if(venturi_bounce(&p)){
                wall_lut[walli].x = x;
                wall_lut[walli].y = y;
                walli++;
            }else{
                // Count particle if in bounds
                int s = get_section(&p);
                if(s >= 0)  section_null_cnts[s] ++;
            }
        }
    }
    // Iterate through the sections to compute areas
    for(int i=0; i<SECTIONS; i++){
        areas[i] = (float)section_null_cnts[i] / (RENDER_WIDTH*HEIGHT/SECTIONS);
        printf("Section %d area:%.3f\n", i, areas[i]);
    }
    wl_actual = walli;
    // Clear the momentum sums
    memset(momentum_sums, 0, sizeof(momentum_sums));
}

// Return Y intercept of particular x location in renderspace
inline int y_intercept(int x){
    particle p = {.vx=0, .vy=0, .x=x, .y=0};
    while(venturi_bounce(&p))  p.y ++;
    return p.y;
}

#define SECT_W  (RENDER_WIDTH/SECTIONS)
#define SECT_W2 (SECT_W*SECT_W)
// Function to determine bounce area in each section for px momentum computation
void get_bounce_areas(float *bounce_areas){
    for(int s=0; s<SECTIONS; s++){
        float xp = SECT_W*s;
        int lb = floorf(UNREND + xp);
        int rb = floorf((UNREND+SECT_W-1.0) + xp);
        int y0 = y_intercept(lb);
        int y1 = y_intercept(rb);
        float slope = (float)(y1-y0)/SECT_W;
        bounce_areas[s] = sqrtf(SECT_W2*(1.0+slope*slope))/SECT_W;
        printf("Sect %d bounce area %.3f\n", s, bounce_areas[s]);
    }
    memset(momentum_sums, 0, sizeof(momentum_sums));
}

// Function to actually render the saved points
void render_walls(SDL_Renderer* renderer, const uint8_t r, const uint8_t g, const uint8_t b){
    // Set render color
    SDL_SetRenderDrawColor(renderer, r, g, b, 255);
    // Render all the points
    SDL_RenderDrawPoints(renderer, wall_lut, wl_actual);
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
static float max_vx = 0.05;
void render_particles_xspd(SDL_Renderer* renderer, particle* particles) {
    max_vx *= 0.99;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        int s = get_section(particles+i);
        float spd = 0;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS) && s>=0){
            // Compute particle speed color index and clamp index
            float vx = particles[i].vx - vx_means[s];
            float vy = particles[i].vy;
            spd = sqrtf(vx*vx + vy*vy);
            int ci = (int)(255.0*(1.0 - spd/max_vx));
            if(ci<0)  ci=0;
            else if(ci>255)  ci=255;
            // Colorize based on particle speed
            SDL_SetRenderDrawColor(renderer, color_lut[ci][0], color_lut[ci][1], color_lut[ci][2], 255);
            SDL_Point pts[circle_pixcnt];
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
            // Find max/min speed
            max_vx = fmaxf(max_vx, spd);
        }
    }
}

// Y speed particle rendering
static float max_vy = 1;
void render_particles_yspd(SDL_Renderer* renderer, particle* particles) {
    float max = 0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int px = (int)particles[i].x - UNREND;
        int py = (int)particles[i].y;
        if(px>RADIUS && px<=(RENDER_WIDTH-RADIUS)){
            float cif = particles[i].vy/max_vy;
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
        // Save max vy
        if(fabsf(particles[i].vy) > max_vy)  max_vy = particles[i].vy;
    }
    max_vy = max;
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

// Function to render particles using a section coloration table (normalized, scaled 0-1.0)
void render_particles_section(SDL_Renderer* renderer, particle* particles, float *section_table_normalized) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int pind = get_section(particles+i);
        if(pind >= 0){
            int ci = (int)floorf(255.99*(1.0-section_table_normalized[pind]));
            if(ci<0)  ci = 0;
            else if(ci>255)  ci=255;
            if(i==0){
                // Draw one particle in magenta always
                SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
            }else{
                // Otherwise sraw by color LUT
                SDL_SetRenderDrawColor(renderer, color_lut[ci][0], color_lut[ci][1], color_lut[ci][2], 255);
            }
            SDL_Point pts[circle_pixcnt];
            float px = particles[i].x - UNREND;
            float py = particles[i].y;
            for(int c=0; c<circle_pixcnt; c++){
                pts[c].x = px + circle_lut[c].x;
                pts[c].y = py + circle_lut[c].y;
            }
            SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
        }
    }
}

uint64_t getus(){
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    uint64_t tm = (uint64_t)now.tv_sec*1000000 + (uint64_t)now.tv_nsec/1000;
    return tm;
}

SDL_Rect inline get_text_size(SDL_Surface *surf, TTF_Font *font, char *str){
    SDL_Color n = {0};
    surf = TTF_RenderUTF8_Blended(font, str, n);
    SDL_Rect r = {.x=0, .y=0, .w=surf->w, .h=surf->h};
    SDL_FreeSurface(surf);
    return r;
}

void inline render_text(SDL_Renderer* renderer, SDL_Surface *surf, SDL_Texture *tex, TTF_Font *font, SDL_Color color, int x, int y, char *str){
    surf = TTF_RenderUTF8_Blended(font, str, color);
    tex = SDL_CreateTextureFromSurface(renderer, surf);
    int w = surf->w;
    SDL_Rect txt_box = { .x=x - w/2, .y=y, .w=w, .h=surf->h};
    SDL_RenderCopy(renderer, tex, NULL, &txt_box);
    SDL_FreeSurface(surf);
    SDL_DestroyTexture(tex);
}

// Copy rectangle of pixels into another renderer
void chunk_copy(SDL_Renderer* renderer, SDL_Rect in_chunk, SDL_Rect out_chunk){
    SDL_Surface *surf = SDL_CreateRGBSurfaceWithFormat(0, in_chunk.w, in_chunk.h, 32, SDL_PIXELFORMAT_RGBA32);
    SDL_RenderReadPixels(renderer, &in_chunk, SDL_PIXELFORMAT_RGBA32, surf->pixels, surf->pitch);
    SDL_Texture *tex = SDL_CreateTextureFromSurface(renderer, surf);
    SDL_FreeSurface(surf);
    out_chunk.h = in_chunk.h;
    out_chunk.w = in_chunk.w;
    SDL_RenderCopy(renderer, tex, NULL, &out_chunk);
    SDL_DestroyTexture(tex);
}

#define THROAT_LEN  RENDER_WIDTH*(THROAT_POS-DIVERGE_START)
//main function for rendering and updating simulation
int main(int argc, char* argv[]) {

    setpriority(PRIO_PROCESS, 0, -20);

    // Get setup defaults and stuff
    char rc = get_file_info(&dt, &exit_px);

    int wall_view_w = fminf(THROAT_LEN, THROAT_WIDTH);

    //initialize SDL
    SDL_Init(SDL_INIT_VIDEO);
    TTF_Init();
    SDL_Window* window = SDL_CreateWindow("Particle Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, RENDER_WIDTH, HEIGHT+THROAT_WIDTH+TXT_HEIGHT*2, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    
    // Pre-compute pixel locations for each rendered circle
    generate_circle(circle_lut, &circle_pixcnt, (int)RADIUS);
    printf("\n");

    // Find the section containing the throat
    particle p0 = {.vx=0, .vy=0, .x=(int)(UNREND+RENDER_WIDTH*THROAT_POS), .y=(int)(HEIGHT*0.5)};
    int throat_sect = get_section(&p0);

    // Build side walls and compute relative section areas
    build_wall(section_areas);
    float bounce_areas[SECTIONS];
    get_bounce_areas(bounce_areas);

    // Init particles
    particle particles[NUM_PARTICLES];
    init_particles(particles, section_areas);

    SDL_GL_SetSwapInterval(VSYNC_PARAM);


    int pupdate=0;
    int kick=0;

    int64_t now, then;
    int64_t tm0 = getus();
    int64_t tm1 = getus();
    int64_t tm_render_particles=0, tm_render_walls=0, tm_present=0, tm_simulate=0, tm_full = 0;
    int tcnt = 0;

    float looptm = 10.0;

    int pcnt = 0;

    float speed[SECTIONS];
    float temperature[SECTIONS];
    float density[SECTIONS];
    float pressure[SECTIONS];
    float mach[SECTIONS];
    float mpx[SECTIONS];

    float t_speed[SECTIONS];
    float t_temperature[SECTIONS];
    float t_density[SECTIONS];
    float t_pressure[SECTIONS];
    float t_mach[SECTIONS];
    float t_mpx[SECTIONS];
    float avg_speed = 0;

    uint32_t loop_cnt = 0;
    float error_integral = 0;
    float last_error = 0;
    float mass_flow_rate = 0;

    TTF_Font *font = TTF_OpenFont("Vremyafwf-Rp1e.ttf", TXT_HEIGHT);
    TTF_Font *big_font = TTF_OpenFont("Vremyafwf-Rp1e.ttf", TXT_HEIGHT*3/2);
    if(!font)  printf("ERROR: %s\n\n", TTF_GetError());
    SDL_Color white = {255, 255, 255, 255};
    SDL_Color magenta = {255, 0, 255, 255};
    SDL_Color yellow = {255, 255, 0, 255};
    SDL_Surface *surf = TTF_RenderUTF8_Blended(font, "Hello", white);
    SDL_Texture *tex;

    float pmax=0, smax=0, dmax=0, tmax=0, mmax=0, mpmax=0;

    float throat_x = RENDER_WIDTH*THROAT_POS + THROAT_WIDTH*0.5;
    float in_x = RENDER_WIDTH - THROAT_WIDTH;
    float long_bp = 0;

    //main loop for simulation
    while (1) {
        now = getus();
        tm_full += now - tm0;
        tm0 = now;

        
        SDL_Event event;
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)  break;
        // Clear screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        // Render walls in gray
        render_walls(renderer, 60, 60, 60);

        // Render particles based on the readut from the start/stop file
        switch (rc){
            case 't':
                render_particles_section(renderer, particles, temperature);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Temperature)");
                break;
            case 'T':
                render_particles_section(renderer, particles, temperature);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Temperature)");
                break;
            case 'p':
                render_particles_section(renderer, particles, pressure);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Pressure)");
                break;
            case 'P':
                render_particles_section(renderer, particles, pressure);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Pressure)");
                break;
            case 'd':
                render_particles_section(renderer, particles, density);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Density)");
                break;
            case 'D':
                render_particles_section(renderer, particles, density);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Density)");
                break;
            case 's':
                render_particles_xspd(renderer, particles);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Speed)");
                break;
            case 'S':
                render_particles_xspd(renderer, particles);
                render_text(renderer, surf, tex, big_font, white, RENDER_WIDTH*(THROAT_POS+DIVERGE_START)/2, HEIGHT-3*TXT_HEIGHT, "(Color by Speed)");
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
        double tsimd = (double)(now-tm1)/1000.0;
        tm1 = now;

        // Render input chunk window
        SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
        float in_spd = (vx_means[SECTIONS-1] + vx_means[SECTIONS-2])/2.0 * dt * tsimd;
        in_x = (in_x < (RENDER_WIDTH - THROAT_WIDTH*1.5)) ? (RENDER_WIDTH - THROAT_WIDTH) : (in_x + in_spd);
        SDL_Rect in_chunk = {.h=THROAT_WIDTH, .w=THROAT_WIDTH, .y=(HEIGHT-THROAT_WIDTH)/2, .x=roundf(in_x)};
        SDL_RenderDrawRect(renderer, &in_chunk);
        SDL_Rect out_loc;
        out_loc.y = HEIGHT + TXT_HEIGHT*2;  out_loc.x = RENDER_WIDTH-THROAT_WIDTH;
        chunk_copy(renderer, in_chunk, out_loc);
        render_text(renderer, surf, tex, font, white, out_loc.x+THROAT_WIDTH/2, HEIGHT+TXT_HEIGHT/2, "Input Bulk View");
        
        // Render throat chunk window
        SDL_SetRenderDrawColor(renderer, 255, 255, 0, 255);
        float throat_spd = (vx_means[throat_sect-1] + vx_means[throat_sect])/2.0 * dt * tsimd;
        throat_x = (throat_x < (RENDER_WIDTH*THROAT_POS - THROAT_WIDTH*3/4)) ? (RENDER_WIDTH*THROAT_POS - THROAT_WIDTH/4) : (throat_x + throat_spd);
        SDL_Rect throat_chunk = {.h=THROAT_WIDTH, .w=THROAT_WIDTH, .y=(HEIGHT-THROAT_WIDTH)/2, .x=roundf(throat_x)};
        SDL_RenderDrawRect(renderer, &throat_chunk);
        out_loc.x -= THROAT_WIDTH*1.5;
        chunk_copy(renderer, throat_chunk, out_loc);
        render_text(renderer, surf, tex, font, white, out_loc.x+THROAT_WIDTH/2, HEIGHT+TXT_HEIGHT/2, "Throat Bulk View");

        // Render the bounce subsections
        SDL_Rect throat_wall = {.x=RENDER_WIDTH*DIVERGE_START, .y=(HEIGHT-THROAT_WIDTH)/2, .w=wall_view_w, .h=wall_view_w/2};
        out_loc.w = throat_wall.w;  out_loc.h=throat_wall.h;
        out_loc.x -= (THROAT_WIDTH*0.5 + wall_view_w);
        chunk_copy(renderer, throat_wall, out_loc);
        int txt_x = out_loc.x+(get_text_size(surf, font, "(T)").w)*0.7 + wall_view_w;
        render_text(renderer, surf, tex, font, yellow, txt_x, out_loc.y+TXT_HEIGHT/4, "(T)");
        SDL_SetRenderDrawColor(renderer, yellow.r, yellow.g, yellow.b, 255);
        SDL_RenderDrawRect(renderer, &out_loc);
        throat_wall.x=RENDER_WIDTH-wall_view_w, throat_wall.y=0, throat_wall.w=wall_view_w, throat_wall.h=wall_view_w/2;
        out_loc.y += wall_view_w/2;
        chunk_copy(renderer, throat_wall, out_loc);
        render_text(renderer, surf, tex, font, white, out_loc.x+wall_view_w/2, HEIGHT+TXT_HEIGHT/2, "Bounce View");
        render_text(renderer, surf, tex, font, magenta, txt_x, out_loc.y+TXT_HEIGHT/4, "(I)");
        SDL_SetRenderDrawColor(renderer, magenta.r, magenta.g, magenta.b, 255);
        SDL_RenderDrawRect(renderer, &out_loc);

        // Render status text on main display
        for(int s=0; s<SECTIONS; s+=1){
            float xloc = (s+0.5)*SECT_W;
            char mstr[10], pstr[10], tstr[10];
            sprintf(mstr, "M%.2f", mach[s]);
            SDL_Rect mt = get_text_size(surf, font, mstr);
            sprintf(tstr, "T%.2f", t_temperature[s]/tmax);
            SDL_Rect tt = get_text_size(surf, font, tstr);
            sprintf(pstr, "P%.2f", t_pressure[s]/pmax);
            SDL_Rect pt = get_text_size(surf, font, pstr);
            int maxw = (mt.w > pt.w) ? mt.w : pt.w;
            maxw = (tt.w > maxw) ? tt.w : maxw;
            SDL_Rect bg = { .x=xloc-maxw/2-TXT_HEIGHT/8, .y=15-TXT_HEIGHT/8, .w=maxw+TXT_HEIGHT/4, .h=15+TXT_HEIGHT+pt.h+tt.h};
            SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 175);
            SDL_RenderFillRect(renderer, &bg);
            render_text(renderer, surf, tex, font, white, xloc, 15, mstr);
            render_text(renderer, surf, tex, font, white, xloc, 15+TXT_HEIGHT, tstr);
            render_text(renderer, surf, tex, font, white, xloc, 15+2*TXT_HEIGHT, pstr);
        }

        // Update renderers
        SDL_RenderPresent(renderer);
        
        if(kick >= KICK_PERIOD){
            kick=0;
            printf("KICK:\n");
            fflush(stdout);
        }
        kick++;
        
        if(pupdate>=PFRAMES){
            // Compute average speed
            avg_speed = speed_sum/speed_cnt;
            // Compute mass flow rate
            mass_flow_rate = 0.5*mass_flow_rate + (0.5/PFRAMES)*((float)mass_flow_sum);
            if(isnan(mass_flow_rate))  mass_flow_rate = 0.0;
            speed_sum = 0;
            speed_cnt = 0;
            mass_flow_cnt = 0;
            mass_flow_sum = 0;

            pcnt++;
            rc = get_file_info( &dt, &exit_px);
            pupdate=0;

            float pmin=0, smin=0, dmin=0, tmin=0, mpmin=0;
            pmax=0, smax=0, dmax=0, tmax=0, mmax=0, mpmax=0;

            float min_press = 0;
            int min_press_sec = 0;
            // Adjust each section by scaling factors
            for(int i=0; i<SECTIONS; i++){
                // X speed given per particle
                vx_means[i] = vx_sums[i] / pcnts[i];
                float spd = fabsf(vx_means[i]);
                // Temperature proportional to KE per particle
                float temp = e_sums[i]/pcnts[i];
                // Mach number is ratio of spd to spd of sound (2D case, c = sqrt(2/3) * Vrms)
                t_mach[i] = 1.225 * spd / sqrtf(temp);
                // Smooth particle avg speed
                t_speed[i] = (1.0-SMOOTHY)*t_speed[i] + SMOOTHY*spd;
                // Smooth particle temperature
                t_temperature[i] = (1.0-SMOOTHY)*t_temperature[i] + SMOOTHY*temp;
                // Density given per unit area
                float dense = pcnts[i]/section_areas[i];
                t_density[i] = (1.0-SMOOTHY)*t_density[i] + SMOOTHY*dense;
                // Ideal gas law pressure is proportional to temp*density
                t_pressure[i] = (1.0-SMOOTHY)*t_pressure[i] + SMOOTHY*temp*dense;
                // Momentum pressure
                t_mpx[i] = (1.0-SMOOTHY)*t_mpx[i] + SMOOTHY*(momentum_sums[i]/bounce_areas[i]);
                // Clear summations for next call
                vx_sums[i] = 0;
                e_sums[i] = 0;
                pcnts[i] = 0;
                momentum_sums[i] = 0;
                // Get minimums
                if(i==0){
                    pmin=t_pressure[0];
                    smin=t_speed[0];
                    dmin=t_density[0];
                    tmin=t_temperature[0];
                    mpmin=t_mpx[0];
                }else{
                    pmin = fminf(t_pressure[i], pmin);
                    smin = fminf(t_speed[i], smin);
                    dmin = fminf(t_density[i], dmin);
                    tmin = fminf(t_temperature[i], tmin);
                    mpmin = fminf(t_mpx[i], mpmin);
                }
                // Compute the maximums for each section
                pmax = fmaxf(t_pressure[i], pmax);
                smax = fmaxf(t_speed[i], smax);
                dmax = fmaxf(t_density[i], dmax);
                tmax = fmaxf(t_temperature[i], tmax);
                mpmax = fmaxf(t_mpx[i], mpmax);
            }

            // Rescale sections to relative values
            int s_mm = 0;
            for(int i=0; i<SECTIONS; i++){
                // Compute relative output metrics
                pressure[i] =   (t_pressure[i]-pmin)/(pmax-pmin);
                speed[i] =      (t_speed[i]-smin)/(smax-smin);
                density[i] =    (t_density[i]-dmin)/(dmax-dmin);
                temperature[i]= (t_temperature[i]-tmin)/(tmax-tmin);
                mpx[i] =        (t_mpx[i]-mpmin)/(mpmax-mpmin);
                // Smooth Mach output
                mach[i] = (1.0-SMOOTHY)*mach[i] + SMOOTHY*t_mach[i];
                // Get raw max Mach
                if(t_mach[i] > mmax){
                    mmax = t_mach[i];
                    s_mm = i;
                }
            }

            printf("SPEEDS: %d", (int)(speed[0]*1000));
            for(int i=1; i<SECTIONS; i++)  printf(" %d", (int)(speed[i]*1000));
            printf("\n");
            printf("PRESSURES: %d", (int)(mpx[0]*1000));
            for(int i=1; i<SECTIONS; i++)  printf(" %d", (int)(mpx[i]*1000));
            printf("\n");
            printf("DENSITIES: %d", (int)(density[0]*1000));
            for(int i=1; i<SECTIONS; i++)  printf(" %d", (int)(density[i]*1000));
            printf("\n");
            printf("TEMPS: %d", (int)(temperature[0]*1000));
            for(int i=1; i<SECTIONS; i++)  printf(" %d", (int)(temperature[i]*1000));
            printf("\n");
            printf("IGLP: %d", (int)(pressure[0]*1000));
            for(int i=1; i<SECTIONS; i++)  printf(" %d", (int)(pressure[i]*1000));
            printf("\n");
            printf("AREA: 0");
            for(int i=1; i<SECTIONS; i++)  printf(" 0");
            printf("\n");
            printf("MACH: %d", (int)(mach[0]*1000));
            for(int i=1; i<SECTIONS; i++)  printf(" %d", (int)(mach[i]*1000));
            printf("\n");
            printf("ENTROPY: 0");
            for(int i=1; i<SECTIONS; i++)  printf(" 0");
            printf("\n");
            looptm = (float)tm_full / 1000.0 / tcnt;
            float simtm = (float)tm_simulate / 1000.0 / tcnt;
            float rentm = ((float)tm_render_particles + (float)tm_render_walls) / 1000.0 / tcnt;
            float prestm = (float)tm_present / 1000.0 / tcnt;
            tcnt=0; tm_full=0; tm_render_particles=0; tm_render_walls=0; tm_simulate=0; tm_present=0;
            int null_cnt = 0;
            int nan_cnt = 0;
            double xum = 0;
            for(int i=0; i<NUM_PARTICLES; i++){
                float x = particles[i].x;
                null_cnt += (x < 0);
                xum += x;
                nan_cnt += isnan(x);
            }
            xum = xum / NUM_PARTICLES;
            printf("SLOWDOWN:%c entry%d/%d/%.2f exit%d/%d/%.2f, null%d, FPS:%.1f, maxM:%.6f throat%.4f in%.4f\n", rc, entry_cnt, init_cnt, entry_inject_temp, exit_cnt, out_cnt_targ, exit_inject_temp, null_cnt, 1000.0/looptm, mmax, throat_spd, in_spd);
            fflush(stdout);
        }
        pupdate++;

        // Simulate for elapsed time from last update call
        simulate_timestep(particles, tsimd);
        tcnt++;
    }

    //clean up SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}