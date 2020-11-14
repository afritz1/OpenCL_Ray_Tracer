#ifndef MAIN_H
#define MAIN_H

#include <SDL\SDL.h>
#include <CL\cl.h>

#define FALSE 0
#define TRUE !FALSE

#define IS_RESIZABLE FALSE
#define IS_OPENGL FALSE
#define IS_FULLSCREEN FALSE
#define IS_HARDWARE_SURFACE TRUE
#define IS_DOUBLE_BUFFERED TRUE
#define SCREEN_FLAGS \
	((IS_RESIZABLE ? SDL_RESIZABLE : 0) | \
	(IS_OPENGL ? SDL_OPENGL : 0) | \
	(IS_FULLSCREEN ? SDL_FULLSCREEN : 0) | \
	(IS_HARDWARE_SURFACE ? SDL_HWSURFACE : SDL_SWSURFACE) | \
	(IS_DOUBLE_BUFFERED ? SDL_DOUBLEBUF : 0)) 

#define SCREEN_TITLE "OpenCL C Ray Tracing! WASD/Space/LShift to move, F/V to change speed, mouse to look"
#define SCREEN_WIDTH ((short)1280)
#define SCREEN_HEIGHT ((short)800)
#define SCREEN_ASPECT ((float)(SCREEN_WIDTH) / (float)(SCREEN_HEIGHT))
#define SCREEN_FPS ((char)60)
#define SCREEN_DELAY (1000 / SCREEN_FPS)
#define SCREEN_BPP ((char)32)

#define MAX_PIXEL_SIZE 1
#define SUPER_SAMPLES 1

#define WORLD_SIZE 12
#define FOG_DIST 2000
#define MAX_REC_DEPTH 30 /* Probably excessive, but it reduces the weird end-of-recursion flat colors. */

#define SPHERE_COUNT 10
#define CUBOID_COUNT 10

#define POINT_LIGHT_COUNT 0
#define DISTANT_LIGHT_COUNT 1
#define AREA_LIGHT_COUNT 0
#define LIGHT_SAMPLES 1 /* For area Light */

#define RENDER_WIDTH (SCREEN_WIDTH / MAX_PIXEL_SIZE)
#define RENDER_HEIGHT (SCREEN_HEIGHT / MAX_PIXEL_SIZE)
#define RENDER_ASPECT ((float)(RENDER_WIDTH) / (float)(RENDER_HEIGHT))

/* Max allowed shapes and lights */
#define MAX_SHAPES 10000
#define MAX_LIGHTS 50

static char doneRendering = FALSE;

// -----------------------
// Vector3 typedef
// -----------------------

typedef struct Vector3
{
	float x, y, z;
} Vector3;

// -----------------------
// Ray typedef
// -----------------------

#define INIT_RAY_DEPTH 1

typedef struct Ray
{
	struct Vector3 point, direction;
	unsigned char depth;
} Ray;

// -----------------------
// Camera typedef
// -----------------------

typedef struct Camera
{
	struct Vector3 eye, forward, right, up;
	float aspect;
	unsigned char superSamples;
} Camera;

// -----------------------
// Light typedef
// -> A light will be queried for its id, and that
//    will determine what to do with its members.
// -----------------------

#define POINT_LIGHT_ID 1
#define DISTANT_LIGHT_ID 2
#define AREA_LIGHT_ID 3

typedef struct Light
{
	unsigned char id;
	struct Vector3 color;
	unsigned short samples;

	/* PointLight, AreaLight */
	struct Vector3 point;

	/* DistantLight */
	struct Vector3 direction;

	/* AreaLight */
	float radius;
} Light;

// -----------------------
// Material typedef
// -> A material will be queried for its id, and that
//    will determine what to do with its members.
// -----------------------

#define FLAT_ID 1
#define PHONG_ID 2

typedef struct Material
{
	unsigned char id;

	/* Flat, Phong */
	struct Vector3 color;

	/* Phong */
	struct Vector3 specularColor;
	float ambient, diffuse, specular;
	unsigned short shiny;
	unsigned char shadows, reflective;
} Material;

// -----------------------
// Shape typedef
// -> A shape will be queried for its id, and that
//    will determine what to do with its members.
// -----------------------

#define SPHERE_ID 1
#define CUBOID_ID 2

typedef struct Shape
{
	unsigned char id;
	struct Vector3 point;
	struct Material material;

	/* Sphere */
	float radius;

	/* Cuboid */
	float width, height, depth;
} Shape;

typedef struct Intersection
{
	float t;
	struct Vector3 point, normal;
	const struct Shape *shape;
} Intersection;

// -----------------------
// World typedef
// -----------------------

typedef struct World
{
	struct Camera camera;
	struct Vector3 backColor, fogColor;
	unsigned short size, maxRecDepth;
	float fogDist;

	struct Shape shapes[MAX_SHAPES];
	unsigned short shapeCount;

	struct Light lights[MAX_LIGHTS];
	unsigned short lightCount;
} World;

// -----------------------
// OpenCL declarations
// -----------------------

typedef struct OpenCL_RT
{
	char *kernelFilename;
	char *kernelString;

	cl_int status;
	cl_platform_id platform;
	cl_device_id device;
	cl_context context;
	cl_command_queue commandQueue;
	cl_program program;
	cl_kernel kernel;

	cl_mem toDoFrame;
	cl_mem worldBuffer;
	cl_uint doneFrame[RENDER_WIDTH * RENDER_HEIGHT];
	size_t globalWorkSize[2];

	size_t startTime;
	size_t endTime;

} OpenCL_RT;

#endif