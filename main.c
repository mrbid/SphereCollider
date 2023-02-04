/*
    James William Fletcher (github.com/mrbid)
        Feb 2023

    Keys:

        O = Orbit in.
        N = New simulation.
        F = FPS to console.
        
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/file.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/time.h>

#pragma GCC diagnostic ignored "-Wunused-result"

#define uint GLushort
#define sint GLshort
#define f32 GLfloat

#include "inc/gl.h"
#define GLFW_INCLUDE_NONE
#include "inc/glfw3.h"

#ifndef __x86_64__
    #define NOSSE
#endif

// uncommenting this define will enable the MMX random when using fRandFloat (it's a marginally slower)
#define SEIR_RAND

#include "inc/esAux2.h"

#include "inc/res.h"
#include "inc/low.h"

//*************************************
// globals
//*************************************
GLFWwindow* window;
uint winw = 1024;
uint winh = 768;
double t = 0;   // time
f32 dt = 0;     // delta time
double fc = 0;  // frame count
double lfct = 0;// last frame count time
double lc = 0;  // logic count
double llct = 0;// last logic count time
f32 aspect;
double x,y,lx,ly;
double rww, ww, rwh, wh, ww2, wh2;
double uw, uh, uw2, uh2; // normalised pixel dpi

// render state id's
GLint projection_id;
GLint modelview_id;
GLint position_id;
GLint lightpos_id;
GLint solidcolor_id;
GLint color_id;
GLint opacity_id;

// render state matrices
mat projection;
mat view;
mat model;
mat modelview;

// render state inputs
vec lightpos = {0.f, 0.f, 0.f};

// models
ESModel mdlSphere;

// game vars
#define FAR_DISTANCE 1000.f
uint RENDER_PASS = 0;
double st=0; // start time
char tts[32];// time taken string
f32 cd = 0.f;
f32 cd2 = 0.f;

// camera vars
uint focus_cursor = 1;
double sens = 0.001f;
f32 xrot = 2.6f;
f32 yrot = 1.f;
f32 zoom = -2.33f;

// sim
#define MAX_SPHERES 1024
f32 SPHERE_SCALE = 0.02f;
f32 SPHERE_SPEED = 0.001f;

typedef struct{vec pos, dir; unsigned char c;} sphere;
sphere spheres[MAX_SPHERES];

//*************************************
// utility functions
//*************************************
void timestamp(char* ts)
{
    const time_t tt = time(0);
    strftime(ts, 16, "%H:%M:%S", localtime(&tt));
}
void timeTaken(uint ss)
{
    if(ss == 1)
    {
        const double tt = t-st;
        if(tt < 60.0)
            sprintf(tts, "%.0f Sec", tt);
        else if(tt < 3600.0)
            sprintf(tts, "%.2f Min", tt * 0.016666667);
        else if(tt < 216000.0)
            sprintf(tts, "%.2f Hr", tt * 0.000277778);
        else if(tt < 12960000.0)
            sprintf(tts, "%.2f Days", tt * 0.00000463);
    }
    else
    {
        const double tt = t-st;
        if(tt < 60.0)
            sprintf(tts, "%.0f Seconds", tt);
        else if(tt < 3600.0)
            sprintf(tts, "%.2f Minutes", tt * 0.016666667);
        else if(tt < 216000.0)
            sprintf(tts, "%.2f Hours", tt * 0.000277778);
        else if(tt < 12960000.0)
            sprintf(tts, "%.2f Days", tt * 0.00000463);
    }
}
float urandf()
{
    static const float RECIP_FLOAT_UINT64_MAX = 1.f/(float)UINT64_MAX;
    int f = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    uint64_t s = 0;
    read(f, &s, sizeof(uint64_t));
    close(f);
    return ((float)s) * RECIP_FLOAT_UINT64_MAX;
}
float urandfc()
{
    static const float RECIP_FLOAT_UINT64_MAX = 2.f/(float)UINT64_MAX;
    int f = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    uint64_t s = 0;
    read(f, &s, sizeof(uint64_t));
    close(f);
    return (((float)s) * RECIP_FLOAT_UINT64_MAX)-1.f;
}
uint64_t urand64()
{
    int f = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    uint64_t s = 0;
    ssize_t result = read(f, &s, sizeof(uint64_t));
    close(f);
    return s;
}
int urand_int()
{
    int f = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    int s = 0;
    read(f, &s, sizeof(int));
    close(f);
    return s;
}
float fRandFloat(const float min, const float max)
{
    return min + randf() * (max-min); 
}
float uRandFloat(const float min, const float max)
{
    return urandf() * (max-min) + min;
}
uint64_t microtime()
{
    struct timeval tv;
    struct timezone tz;
    memset(&tz, 0, sizeof(struct timezone));
    gettimeofday(&tv, &tz);
    return 1000000 * tv.tv_sec + tv.tv_usec;
}

//*************************************
// render functions
//*************************************

void rSphere(f32 x, f32 y, f32 z)
{
    mIdent(&model);
    mTranslate(&model, x, y, z);
    mScale(&model, SPHERE_SCALE, SPHERE_SCALE, SPHERE_SCALE);
    mMul(&modelview, &model, &view);

    glUniformMatrix4fv(modelview_id, 1, GL_FALSE, (f32*) &modelview.m[0][0]);
    glDrawElements(GL_TRIANGLES, low_numind, GL_UNSIGNED_SHORT, 0);
}

//*************************************
// game functions
//*************************************

void newSim()
{
    const int seed = urand_int();
    srand(seed);
    srandf(seed);

    st = t;

    char strts[16];
    timestamp(&strts[0]);
    printf("[%s] Sim Start [%u].\n", strts, seed);
    
    glfwSetWindowTitle(window, "Sphere Collider");
    
    for(uint i = 0; i < MAX_SPHERES; i++)
    {
        vRuvTA(&spheres[i].pos); // random point on inside of unit sphere
        vRuvBT(&spheres[i].dir); // random point on outside of unit sphere
        vNorm(&spheres[i].dir); // hmm or not?
        spheres[i].c = 0.f;
    }
}

//*************************************
// update & render
//*************************************
void main_loop()
{
//*************************************
// camera
//*************************************

    if(focus_cursor == 1)
    {
        glfwGetCursorPos(window, &x, &y);

        xrot += (ww2-x)*sens;
        yrot += (wh2-y)*sens;

        if(yrot > 1.5f)
            yrot = 1.5f;
        if(yrot < 0.5f)
            yrot = 0.5f;

        glfwSetCursorPos(window, ww2, wh2);
    }

    mIdent(&view);
    mTranslate(&view, 0.f, 0.f, zoom);
    mRotate(&view, yrot, 1.f, 0.f, 0.f);
    mRotate(&view, xrot, 0.f, 0.f, 1.f);

//*************************************
// render
//*************************************
    if(RENDER_PASS == 1)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    for(uint i = 0; i < MAX_SPHERES; i++)
    {
        if(RENDER_PASS == 1)
        {
            if(spheres[i].c == 0)
                glUniform3f(color_id, 0.f, 0.f, 1.f);
            else
                glUniform3f(color_id, 1.f, 0.f, 0.f);
        }

        vec inc;
        vMulS(&inc, spheres[i].dir, SPHERE_SPEED);
        vAdd(&spheres[i].pos, spheres[i].pos, inc);

        const f32 mod = vMod(spheres[i].pos);
        if(mod > 1.f)
        {
            if(RENDER_PASS == 1 && mod > 1.01f){glUniform3f(color_id, 1.f, 1.f, 0.f);}

            vec sd = spheres[i].pos;
            vNorm(&sd);

            vReflect(&spheres[i].dir, spheres[i].dir, sd);
            vNorm(&spheres[i].dir); // hmm or not?

            vec ob = spheres[i].pos;
            vNorm(&ob);
            vInv(&ob);

            vec inc;
            if(mod < 1.03f)
                vMulS(&inc, spheres[i].dir, (mod-1.f)+SPHERE_SPEED);
            else
                vMulS(&inc, ob, (mod-1.f)+SPHERE_SPEED);
            vAdd(&spheres[i].pos, spheres[i].pos, inc);
        }

        for(uint j = 0; j < MAX_SPHERES; j++)
        {
            if(j == i){continue;} // dont collide with self

            f32 d = vDistSq(spheres[i].pos, spheres[j].pos);
            if(d < cd2)
            {
                if(RENDER_PASS == 1){glUniform3f(color_id, 1.f, 0.f, 0.f);}

                // vec sdj = spheres[j].pos;
                // vNorm(&sdj);
                // vReflect(&spheres[i].dir, sdj, spheres[i].dir);
                // vNorm(&spheres[i].dir);
                
                // reflect the ball direction
                vReflect(&spheres[i].dir, spheres[j].dir, spheres[i].dir);
                vNorm(&spheres[i].dir);

                // increment the ball to a non-inersecting distance in the new direction
                vec inc;
                d = vDist(spheres[i].pos, spheres[j].pos);
                vMulS(&inc, spheres[i].dir, (cd-d)+SPHERE_SPEED);
                vAdd(&spheres[i].pos, spheres[i].pos, inc);

                // char strts[16];
                // timestamp(&strts[0]);
                // printf("[%s] Collision %f.\n", strts, d);
            }
        }
        
        if(RENDER_PASS == 1)
            rSphere(spheres[i].pos.x, spheres[i].pos.y, spheres[i].pos.z);
    }

    if(RENDER_PASS == 1)
        glfwSwapBuffers(window);
}

//*************************************
// Input Handelling
//*************************************
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // control
    if(action == GLFW_PRESS)
    {
        // new
        if(key == GLFW_KEY_N)
        {
            // end
            timeTaken(0);
            char strts[16];
            timestamp(&strts[0]);
            printf("[%s] Sim End.\n", strts);
            printf("[%s] Time-Taken: %s or %g Seconds\n\n", strts, tts, t-st);
            
            // new
            newSim();
        }

        // show average fps
        else if(key == GLFW_KEY_F)
        {
            if(t-lfct > 2.0)
            {
                char strts[16];
                timestamp(&strts[0]);
                printf("[%s] FPS: %g\n", strts, fc/(t-lfct));
                printf("[%s] LPS: %g\n", strts, lc/(t-llct));
                lfct = t;
                fc = 0;
                llct = t;
                lc = 0;
            }
        }

        // orbit in
        else if(key == GLFW_KEY_O)
        {
            for(uint i = 0; i < MAX_SPHERES; i++)
            {
                vRuvBT(&spheres[i].pos); // random point on outside of unit sphere
                vMulS(&spheres[i].pos, spheres[i].pos, 1.01f+(randf()*0.02f)); // project it away from the unit sphere a little
                vRuvBT(&spheres[i].dir); // random point on outside of unit sphere
                vNorm(&spheres[i].dir); // hmm or not?
                spheres[i].c = 0.f;
            }
        }
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    if(yoffset == -1)
        zoom += 0.06f * zoom;
    else if(yoffset == 1)
        zoom -= 0.06f * zoom;
    
    if(zoom > -2.33f){zoom = -2.33f;}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        if(button == GLFW_MOUSE_BUTTON_LEFT)
        {
            focus_cursor = 1 - focus_cursor;
            if(focus_cursor == 0)
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            else
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
            glfwSetCursorPos(window, ww2, wh2);
            glfwGetCursorPos(window, &ww2, &wh2);
        }
        else if(button == GLFW_MOUSE_BUTTON_4 || button == GLFW_MOUSE_BUTTON_RIGHT)
        {
            if(zoom != -0.3f)
                zoom = -0.3f;
            else
                zoom = -2.33f;
        }
    }
}

void window_size_callback(GLFWwindow* window, int width, int height)
{
    winw = width;
    winh = height;

    glViewport(0, 0, winw, winh);
    ww = winw;
    wh = winh;
    aspect = ww / wh;
    rww = 1/ww;
    rwh = 1/wh;
    ww2 = ww/2;
    wh2 = wh/2;
    uw = (double)aspect / ww;
    uh = 1 / wh;
    uw2 = (double)aspect / ww2;
    uh2 = 1 / wh2;

    mIdent(&projection);
    mPerspective(&projection, 60.0f, aspect, 0.1f, FAR_DISTANCE);
    glUniformMatrix4fv(projection_id, 1, GL_FALSE, (f32*) &projection.m[0][0]);
}

//*************************************
// Process Entry Point
//*************************************
int main(int argc, char** argv)
{
    // allow custom msaa level
    int msaa = 16;
    if(argc >= 2){msaa = atoi(argv[1]);}

    // allow framerate cap
    double maxfps = 144.0;
    if(argc >= 3){maxfps = atof(argv[2]);}

    // help
    printf("----\n");
    printf("Sphere Collider\n");
    printf("----\n");
    printf("James William Fletcher (github.com/mrbid)\n");
    printf("----\n");
    printf("Argv(2): msaa, maxfps\n");
    printf("e.g; ./uc 16 60\n");
    printf("----\n");

    // init glfw
    if(!glfwInit()){printf("glfwInit() failed.\n"); exit(EXIT_FAILURE);}
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_SAMPLES, msaa);
    window = glfwCreateWindow(winw, winh, "Sphere Collider", NULL, NULL);
    if(!window)
    {
        printf("glfwCreateWindow() failed.\n");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    const GLFWvidmode* desktop = glfwGetVideoMode(glfwGetPrimaryMonitor());
    glfwSetWindowPos(window, (desktop->width/2)-(winw/2), (desktop->height/2)-(winh/2)); // center window on desktop
    glfwSetWindowSizeCallback(window, window_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwMakeContextCurrent(window);
    gladLoadGL(glfwGetProcAddress);
    glfwSwapInterval(0); // 0 for immediate updates, 1 for updates synchronized with the vertical retrace, -1 for adaptive vsync

    // set icon
    glfwSetWindowIcon(window, 1, &(GLFWimage){16, 16, (unsigned char*)&icon_image.pixel_data});

    // hide cursor
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);


//*************************************
// projection
//*************************************

    window_size_callback(window, winw, winh);

//*************************************
// bind vertex and index buffers
//*************************************

    esBind(GL_ARRAY_BUFFER, &mdlSphere.vid, low_vertices, sizeof(low_vertices), GL_STATIC_DRAW);
    esBind(GL_ELEMENT_ARRAY_BUFFER, &mdlSphere.iid, low_indices, sizeof(low_indices), GL_STATIC_DRAW);

//*************************************
// compile & link shader programs
//*************************************

    makeLambert();

//*************************************
// configure render options
//*************************************
    
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.13, 0.13, 0.13, 0.0);

    shadeLambert(&position_id, &projection_id, &modelview_id, &lightpos_id, &color_id, &opacity_id);
    glUniformMatrix4fv(projection_id, 1, GL_FALSE, (f32*) &projection.m[0][0]);
    glUniform3f(lightpos_id, lightpos.x, lightpos.y, lightpos.z);
    glUniform1f(opacity_id, 1.0f);

    glBindBuffer(GL_ARRAY_BUFFER, mdlSphere.vid);
    glVertexAttribPointer(position_id, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_id);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mdlSphere.iid);

//*************************************
// execute update / render loop
//*************************************

    // init
    cd = SPHERE_SCALE*1.8f;
    cd2 = cd*cd;
    const double maxlps = 144.0;
    t = glfwGetTime();
    lfct = t;
    dt = 1.0 / (float)maxlps; // fixed timestep delta-time
    newSim();

    // just for launch
    for(uint i = 0; i < MAX_SPHERES; i++)
    {
        vRuvBT(&spheres[i].pos); // random point on outside of unit sphere
        vMulS(&spheres[i].pos, spheres[i].pos, 1.01f+(randf()*0.02f)); // project it away from the unit sphere a little
        vRuvBT(&spheres[i].dir); // random point on outside of unit sphere
        vNorm(&spheres[i].dir); // hmm or not?
        spheres[i].c = 0.f;
    }
    
    // lps accurate event loop
    const double fps_limit = 1.0 / maxfps;
    double rlim = 0.0;
    const useconds_t wait_interval = 1000000 / maxlps; // fixed timestep
    useconds_t wait = wait_interval;
    while(!glfwWindowShouldClose(window))
    {
        usleep(wait);
        t = glfwGetTime();
        glfwPollEvents();

        if(maxfps < maxlps)
        {
            if(t > rlim)
            {
                RENDER_PASS = 1;
                rlim = t + fps_limit; // should be doing this after main_loop() at the very least but it's not a big deal.
                fc++;
            }
            else
                RENDER_PASS = 0;
        }
        else
        {
            RENDER_PASS = 1;
            fc++;
        }

        main_loop();

        // accurate lps
        wait = wait_interval - (useconds_t)((glfwGetTime() - t) * 1000000.0);
        if(wait > wait_interval)
            wait = wait_interval;
        lc++;
    }

    // end
    timeTaken(0);
    char strts[16];
    timestamp(&strts[0]);
    printf("[%s] Sim End.\n", strts);
    printf("[%s] Time-Taken: %s or %g Seconds\n\n", strts, tts, t-st);

    // done
    glfwDestroyWindow(window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
    return 0;
}
