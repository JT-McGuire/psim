#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL2/SDL.h>

#define WIDTH 2000.0 //width of bounding box
#define HEIGHT 1200.0 //height of bounding box
#define WALL_RADIUS 1000.0
#define NECK_HEIGHT 200.0
#define NUM_PARTICLES 1000 //number of particles in the simulation
#define RADIUS 10.0 //radius of each particle
#define DT 0.01 //timestep for simulation

#define INITIAL_TEMP 500.0


#define REPOS 1.0001
#define FORCE 5.0

#define INCOMPRESSIBLE


//structure for a single particle
typedef struct {
    float x; //x position of particle
    float y; //y position of particle
    float vx; //x velocity of particle
    float vy; //y velocity of particle
} particle;

// Returns a random float between 0 and 1
inline float randf(){ return ((float)rand() / RAND_MAX); }

#define BIGCIRC_X (WIDTH/2)
#define BIGCIRC_YH (HEIGHT/2+NECK_HEIGHT/2+WALL_RADIUS)
#define BIGCIRC_YL (HEIGHT/2-NECK_HEIGHT/2-WALL_RADIUS)
#define PUSH_RAD ((WALL_RADIUS+RADIUS)*REPOS)
uint8_t venturi_bounce(particle *p){
    // Check X rectangle bounds
    if(p->x <= RADIUS){
        //p->vx *= -1;
        //p->x = RADIUS*REPOS;
        p->x = WIDTH-RADIUS*REPOS;
        return 1;
    }
    if(p->x >= (WIDTH-RADIUS)){
        p->vx *= -1;
        p->x = WIDTH-RADIUS*REPOS;
        return 1;
    }
    // Check Y rectangle bounds
    if(p->y <= RADIUS){
        p->vy *= -1;
        p->y = RADIUS*REPOS;
        return 1;
    }
    if(p->y >= (HEIGHT-RADIUS)){
        p->vy *= -1;
        p->y = HEIGHT-RADIUS*REPOS;
        return 1;
    }
    // Compute distance from upper part of constriction bound
    float dx = p->x - BIGCIRC_X;
    float dy = p->y - BIGCIRC_YH;
    float dx2 = dx*dx;
    float dist = sqrt(dx2 + dy*dy);
    float nx, ny, s;
    if( dist<(WALL_RADIUS+RADIUS) ){
        // Perform householder reflection of velocity
        nx = dx/dist;
        ny = dy/dist;
        s = nx*p->vx + ny*p->vy;
        p->vx -= nx*2*s;
        p->vy -= ny*2*s;
        //printf("Norm:(%f,%f)", nx, ny);
        //printf(" Orig:(%f,%f)", p->x, p->y);
        p->x = BIGCIRC_X + PUSH_RAD*nx;
        p->y = BIGCIRC_YH + PUSH_RAD*ny;
        //printf(" Bounce:(%f,%f)\n", p->x, p->y);
        return 1;
    }
    dy = p->y - BIGCIRC_YL;
    dist = sqrt(dx2 + dy*dy);
    if( dist<(WALL_RADIUS+RADIUS) ){
        // Perform householder reflection of velocity
        nx = dx/dist;
        ny = dy/dist;
        s = nx*p->vx + ny*p->vy;
        p->vx -= nx*2*s;
        p->vy -= ny*2*s;
        p->x = BIGCIRC_X + PUSH_RAD*nx;
        p->y = BIGCIRC_YL + PUSH_RAD*ny;
        return 1;
    }
    return 0;
}

// Function to check if 
inline uint8_t check_overlap(particle* particles, float x, float y, int occ){
    particle p;
    p.x=x; p.y=y; p.vx=0; p.vy=0;
    if(venturi_bounce(&p)){ return 1; }
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
            particles[i].vx = (randf()*2 - 1)*INITIAL_TEMP;
            particles[i].vy = (randf()*2 - 1)*INITIAL_TEMP;
            i++;
        }
    }
}

// Function to simulate a timestep for each particle
void simulate_timestep(particle* particles, float ideal_dist) {
    // Scan through all particle pair combos
    for (int i = 0; i < NUM_PARTICLES - 1; i++) {
        for (int j = i + 1; j < NUM_PARTICLES; j++) {
            // Compute distance with distance formula
            float dx = particles[i].x - particles[j].x;
            float dy = particles[i].y - particles[j].y;
            float dist = sqrt(dx*dx + dy*dy);
            // Compute normal vector components
            float nx = dx/dist;
            float ny = dy/dist;
            if (dist <= 2*RADIUS) {
                // Compute scalar impulse value
                float k = nx*(particles[j].vx - particles[i].vx) + ny*(particles[j].vy - particles[i].vy);
                // Add the deltaV to each velocity
                particles[i].vx += k*nx;
                particles[i].vy += k*ny;
                particles[j].vx -= k*nx;
                particles[j].vy -= k*ny;
                particles[i].x = particles[j].x + nx*(2*RADIUS*REPOS);
                particles[i].y = particles[j].y + ny*(2*RADIUS*REPOS);
            }
            #ifdef INCOMPRESSIBLE
                // Compute the acceleration from the distance
                float accel = FORCE/(dist*dist);
                float ax = accel*nx;
                float ay = accel*ny;
                particles[i].vx += ax;
                particles[i].vy += ay;
                particles[j].vx -= ax;
                particles[j].vy -= ay;
            #endif
        }
    }

    // Update positions of particles
    for (int i = 0; i < NUM_PARTICLES; i++) {
        particles[i].x += particles[i].vx * DT;
        particles[i].y += particles[i].vy * DT;
    }

    // Check for collisions with bounding box
    for (int i = 0; i < NUM_PARTICLES; i++) {
        venturi_bounce(particles+i);
    }
}

static SDL_Point circle_lut[(int)((RADIUS+1)*(RADIUS+1)*4)];
static int circle_pixcnt=0;

static SDL_Point wall_high[(int)(WALL_RADIUS*WALL_RADIUS*M_PI/2.0)];
static int wall_cnt_h=0;
static SDL_Point wall_low[(int)(WALL_RADIUS*WALL_RADIUS*M_PI/2.0)];
static int wall_cnt_l=0;

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
                if(!use_bounds || (x>=0 && x<WIDTH && y>=0 && y<HEIGHT)){
                    lut[(*pixcnt)].x = x;
                    lut[(*pixcnt)].y = y;
                    (*pixcnt) ++;
                    //printf("Circle point (%d,%d) count %d\n", x, y, circle_pixcnt);
                }
            }
        }
    }
}

//function to render particles on the screen
void render_particles(SDL_Renderer* renderer, particle* particles) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_Point pts[circle_pixcnt];
        for(int c=0; c<circle_pixcnt; c++){
            pts[c].x = (int)particles[i].x + circle_lut[c].x;
            pts[c].y = (int)particles[i].y + circle_lut[c].y;
        }
        SDL_RenderDrawPoints(renderer, pts, circle_pixcnt);
    }
}

void draw_walls(SDL_Renderer* renderer){
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    SDL_RenderDrawPoints(renderer, wall_high, wall_cnt_h);
    SDL_RenderDrawPoints(renderer, wall_low, wall_cnt_l);
}

//main function for rendering and updating simulation
int main(int argc, char* argv[]) {
    //initialize SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("Particle Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    //initialize particles
    particle particles[NUM_PARTICLES];
    init_particles(particles);

    generate_circle(circle_lut, &circle_pixcnt, 0, 0, (int)RADIUS);
    generate_circle(wall_high, &wall_cnt_h, (int)BIGCIRC_X, (int)BIGCIRC_YH, (int)WALL_RADIUS);
    generate_circle(wall_low, &wall_cnt_l, (int)BIGCIRC_X, (int)BIGCIRC_YL, (int)WALL_RADIUS);

    // Compute the mean inter-particle distance assuming even distribution
    float ideal_circ_area = (WIDTH*HEIGHT - wall_cnt_h - wall_cnt_l)*0.9069 / (float)NUM_PARTICLES;
    float part_dist = 2*sqrtf(ideal_circ_area/M_PI);

    printf("Particle Distance Target is %.3f pixels\n\n", part_dist);

    printf("Top big circle center (%f, %f)\n", BIGCIRC_X, BIGCIRC_YH);
    printf("Bottom big circle center (%f, %f)\n", BIGCIRC_X, BIGCIRC_YL);

    //main loop for simulation
    while (1) {
        SDL_Event event;
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
            break;
        }

        //clear screen and render particles
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        render_particles(renderer, particles);
        draw_walls(renderer);
        SDL_RenderPresent(renderer);

        //simulate timestep for particles
        simulate_timestep(particles, part_dist);
        usleep(1000);
    }

    //clean up SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}