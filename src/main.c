#include "main.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// ---------------------------------
// Utility functions
// ---------------------------------

void printVec(const char *name, const Vector3 *v)
{
	printf("%s[%.3f, %.3f, %.3f]\n", name, v->x, v->y, v->z);
}

// ---------------------------------
// HandleInput()
// ---------------------------------

void handleInput(World *world)
{
	int x = 0, y = 0;
	size_t clicked = SDL_GetRelativeMouseState(&x, &y);
	float mouseSpeed = 0.15f*sqrtf((float)(x*x*0.4f + y*y));

#pragma warning(push)
#pragma warning(disable:4800)
	unsigned char *keys = SDL_GetKeyState(NULL);
	char leftClick = clicked & SDL_BUTTON(SDL_BUTTON_LEFT);
	char rightClick = clicked & SDL_BUTTON(SDL_BUTTON_RIGHT);
	char forward = keys[SDLK_w] | keys[SDLK_UP];
	char backward = keys[SDLK_s] | keys[SDLK_DOWN];
	char left = keys[SDLK_a] | keys[SDLK_LEFT];
	char right = keys[SDLK_d] | keys[SDLK_RIGHT];
	char up = keys[SDLK_SPACE];
	char down = keys[SDLK_LSHIFT];
	char turnRight = (x > 0) && (rightClick || leftClick);
	char turnLeft = (x < 0) && (rightClick || leftClick);
	char turnUp = (y < 0) && (rightClick || leftClick);
	char turnDown = (y > 0) && (rightClick || leftClick);
	char faster = keys[SDLK_f];
	char slower = keys[SDLK_v];
	char any = forward | backward | left | right | up | down |
		turnRight | turnLeft | turnUp | turnDown | faster | slower;
#pragma warning(pop)

	Camera camera = world->camera;
	Vector3 globalUp = { 0.0f, 1.0f, 0.0f };
	Vector3 newEye = camera.eye;
	Vector3 newFwd = camera.forward;
	Vector3 newRgt = camera.right;
	Vector3 newUp = camera.up;
	static float velocity = 0.15f;

	if (!any) return;

	if (faster) velocity *= 1.05f;
	if (slower) velocity *= 0.95f;

	newEye.x +=
		(newFwd.x*velocity*forward) -
		(newFwd.x*velocity*backward) +
		(newRgt.x*velocity*right) -
		(newRgt.x*velocity*left) +
		(newUp.x*velocity*up) -
		(newUp.x*velocity*down);
	newEye.y +=
		(newFwd.y*velocity*forward) -
		(newFwd.y*velocity*backward) +
		(newRgt.y*velocity*right) -
		(newRgt.y*velocity*left) +
		(newUp.y*velocity*up) -
		(newUp.y*velocity*down);
	newEye.z +=
		(newFwd.z*velocity*forward) -
		(newFwd.z*velocity*backward) +
		(newRgt.z*velocity*right) -
		(newRgt.z*velocity*left) +
		(newUp.z*velocity*up) -
		(newUp.z*velocity*down);

	newFwd.x +=
		(newRgt.x*(0.05f*mouseSpeed)*turnRight) -
		(newRgt.x*(0.05f*mouseSpeed)*turnLeft) +
		(newUp.x*(0.05f*mouseSpeed)*turnUp) -
		(newUp.x*(0.05f*mouseSpeed)*turnDown);
	newFwd.y +=
		(newRgt.y*(0.05f*mouseSpeed)*turnRight) -
		(newRgt.y*(0.05f*mouseSpeed)*turnLeft) +
		(newUp.y*(0.05f*mouseSpeed)*turnUp) -
		(newUp.y*(0.05f*mouseSpeed)*turnDown);
	newFwd.z +=
		(newRgt.z*(0.05f*mouseSpeed)*turnRight) -
		(newRgt.z*(0.05f*mouseSpeed)*turnLeft) +
		(newUp.z*(0.05f*mouseSpeed)*turnUp) -
		(newUp.z*(0.05f*mouseSpeed)*turnDown);

	float newFwdLengthRecip = 1.0f / sqrtf(
		(newFwd.x*newFwd.x) + (newFwd.y*newFwd.y) + (newFwd.z*newFwd.z));
	newFwd.x *= newFwdLengthRecip;
	newFwd.y *= newFwdLengthRecip;
	newFwd.z *= newFwdLengthRecip;

	newRgt.x = (newFwd.y * globalUp.z) - (globalUp.y * newFwd.z);
	newRgt.y = (globalUp.x * newFwd.z) - (newFwd.x * globalUp.z);
	newRgt.z = (newFwd.x * globalUp.y) - (globalUp.x * newFwd.y);
	float newRgtLengthRecip = 1.0f / sqrtf(
		(newRgt.x*newRgt.x) + (newRgt.y*newRgt.y) + (newRgt.z*newRgt.z));
	newRgt.x *= newRgtLengthRecip * camera.aspect;
	newRgt.y *= newRgtLengthRecip * camera.aspect;
	newRgt.z *= newRgtLengthRecip * camera.aspect;

	newUp.x = (newRgt.y * newFwd.z) - (newFwd.y * newRgt.z);
	newUp.y = (newFwd.x * newRgt.z) - (newRgt.x * newFwd.z);
	newUp.z = (newRgt.x * newFwd.y) - (newFwd.x * newRgt.y);
	float newUpLengthRecip = 1.0f / sqrtf(
		(newUp.x*newUp.x) + (newUp.y*newUp.y) + (newUp.z*newUp.z));
	newUp.x *= newUpLengthRecip;
	newUp.y *= newUpLengthRecip;
	newUp.z *= newUpLengthRecip;

	if (any)
	{
		doneRendering = FALSE;
		camera.eye = newEye;
		camera.forward = newFwd;
		camera.right = newRgt;
		camera.up = newUp;
		world->camera = camera;
	}
}

// ---------------------------------
// Builder functions
// ---------------------------------

Camera lookAt(const Vector3 *eye, const Vector3 *focus, const Vector3 *up,
	float fovY, float aspect)
{
	Camera camera;
	camera.eye = *eye;
	camera.aspect = aspect;
	camera.superSamples = SUPER_SAMPLES;
	float fovYToRad = fovY * (3.1415926536f / 180.0f);
	float imgHalfH = tanf(fovYToRad / 2.0f);
	float imgHalfW = imgHalfH * aspect;

	camera.forward.x = focus->x - eye->x;
	camera.forward.y = focus->y - eye->y;
	camera.forward.z = focus->z - eye->z;

	float forwardLengthRecip = 1.0f / sqrtf(
		(camera.forward.x*camera.forward.x) +
		(camera.forward.y*camera.forward.y) +
		(camera.forward.z*camera.forward.z));

	camera.forward.x *= forwardLengthRecip;
	camera.forward.y *= forwardLengthRecip;
	camera.forward.z *= forwardLengthRecip;

	camera.right.x = (camera.forward.y * up->z) - (up->y * camera.forward.z);
	camera.right.y = (up->x * camera.forward.z) - (camera.forward.x * up->z);
	camera.right.z = (camera.forward.x * up->y) - (up->x * camera.forward.y);

	float rightLengthRecip = 1.0f / sqrtf(
		(camera.right.x*camera.right.x) +
		(camera.right.y*camera.right.y) +
		(camera.right.z*camera.right.z));

	camera.right.x *= rightLengthRecip;
	camera.right.y *= rightLengthRecip;
	camera.right.z *= rightLengthRecip;

	camera.up.x = (camera.right.y * camera.forward.z) - (camera.forward.y * camera.right.z);
	camera.up.y = (camera.forward.x * camera.right.z) - (camera.right.x * camera.forward.z);
	camera.up.z = (camera.right.x * camera.forward.y) - (camera.forward.x * camera.right.y);

	float upLengthRecip = 1.0f / sqrtf(
		(camera.up.x*camera.up.x) +
		(camera.up.y*camera.up.y) +
		(camera.up.z*camera.up.z));

	camera.up.x *= upLengthRecip;
	camera.up.y *= upLengthRecip;
	camera.up.z *= upLengthRecip;

	camera.right.x *= imgHalfW;
	camera.right.y *= imgHalfW;
	camera.right.z *= imgHalfW;
	camera.up.x *= imgHalfH;
	camera.up.y *= imgHalfH;
	camera.up.z *= imgHalfH;

	return camera;
}

Material makeFlat()
{
	Material m;
	m.id = FLAT_ID;
	m.color.x = (float)((rand() % 100) / 100.0f);
	m.color.y = (float)((rand() % 100) / 100.0f);
	m.color.z = (float)((rand() % 100) / 100.0f);
	return m;
}

Material makePhong()
{
	Material m;
	m.id = PHONG_ID;
	m.color.x = (float)((rand() % 100) / 100.0f);
	m.color.y = (float)((rand() % 100) / 100.0f);
	m.color.z = (float)((rand() % 100) / 100.0f);
	m.specularColor.x = (float)((rand() % 100) / 100.0f);
	m.specularColor.y = (float)((rand() % 100) / 100.0f);
	m.specularColor.z = (float)((rand() % 100) / 100.0f);
	m.ambient = 0.25f;
	m.diffuse = 1.00f;
	m.specular = 0.40f;
	m.shiny = 32;
	m.shadows = TRUE;
	m.reflective = TRUE;
	return m;
}

Shape makeSphere(const World *this)
{
	Shape s;
	s.id = SPHERE_ID;
	s.point.x = (float)(rand() % this->size) - (this->size * 0.5f);
	s.point.y = (float)(rand() % this->size) + (this->size * 0.5f);
	s.point.z = (float)(rand() % this->size) - (this->size * 0.5f);
	s.material = makePhong();
	s.radius = (0.5f + (float)((rand() % 100) / 100.0f));
	return s;
}

Shape makeCuboid(const World *this)
{
	Shape s;
	s.id = CUBOID_ID;
	s.point.x = (float)(rand() % this->size) - (this->size * 0.5f);
	s.point.y = (float)(rand() % this->size) + (this->size * 0.5f);
	s.point.z = (float)(rand() % this->size) - (this->size * 0.5f);
	s.material = makePhong();
	s.width = (0.5f + (float)((rand() % 100) / 100.0f));
	s.height = (0.5f + (float)((rand() % 100) / 100.0f));
	s.depth = (0.5f + (float)((rand() % 100) / 100.0f));
	return s;
}

Shape makeCuboidPlane(const World *this)
{
	Shape s = makeCuboid(this);
	s.width = CL_MAXFLOAT;
	s.height = 0.0f;
	s.depth = CL_MAXFLOAT;
	s.point.x = 0.0f;
	s.point.y = 0.0f;
	s.point.z = 0.0f;
	return s;
}

Light makePointLight(const World *this)
{
	Light l;
	l.id = POINT_LIGHT_ID;
	l.samples = 1;
	l.color.x = (float)((rand() % 100) / 100.0f);
	l.color.y = (float)((rand() % 100) / 100.0f);
	l.color.z = (float)((rand() % 100) / 100.0f);
	l.point.x = (float)(rand() % this->size) - (this->size * 0.5f);
	l.point.y = (float)(rand() % this->size) + (this->size * 1.5f);
	l.point.z = (float)(rand() % this->size) - (this->size * 0.5f);
	return l;
}

Light makeDistantLight(const World *this)
{
	Light l;
	l.id = DISTANT_LIGHT_ID;
	l.samples = 1;
	l.color.x = (float)((rand() % 100) / 100.0f);
	l.color.y = (float)((rand() % 100) / 100.0f);
	l.color.z = (float)((rand() % 100) / 100.0f);
	l.direction.x = (float)(rand() % this->size) - (this->size * 0.5f);
	l.direction.y = (float)(rand() % this->size) + (this->size * 0.05f);
	l.direction.z = (float)(rand() % this->size) - (this->size * 0.5f);
	float dirLen = 1.0f / sqrtf(
		(l.direction.x*l.direction.x) +
		(l.direction.y*l.direction.y) +
		(l.direction.z*l.direction.z));
	l.direction.x *= dirLen;
	l.direction.y *= dirLen;
	l.direction.z *= dirLen;
	return l;
}

Light makeAreaLight(const World *this)
{
	Light l;
	l.id = AREA_LIGHT_ID;
	l.samples = LIGHT_SAMPLES;
	l.radius = 0.45f;
	l.color.x = (float)((rand() % 100) / 100.0f);
	l.color.y = (float)((rand() % 100) / 100.0f);
	l.color.z = (float)((rand() % 100) / 100.0f);
	l.point.x = (float)(rand() % this->size) - (this->size * 0.5f);
	l.point.y = (float)(rand() % this->size) + (this->size * 0.5f);
	l.point.z = (float)(rand() % this->size) - (this->size * 0.5f);
	return l;
}

World makeWorld()
{
	World world;

	Vector3 eye = { 8.0f, WORLD_SIZE + 4.0f, 12.0f };
	Vector3 focus = { 0.0f, WORLD_SIZE, 0.0f };
	Vector3 up = { 0.0f, 1.0f, 0.0f };
	world.camera = lookAt(&eye, &focus, &up, 90.0f, SCREEN_ASPECT);

	world.size = WORLD_SIZE;
	world.fogDist = (float)FOG_DIST;
	world.maxRecDepth = MAX_REC_DEPTH;

	world.shapeCount = SPHERE_COUNT + CUBOID_COUNT;
	world.lightCount = POINT_LIGHT_COUNT + DISTANT_LIGHT_COUNT + AREA_LIGHT_COUNT;

	if (world.shapeCount > MAX_SHAPES || world.lightCount > MAX_LIGHTS)
	{
		printf("Error! Too many world objects or lights to construct!\n");
		getchar();
		exit(EXIT_FAILURE);
	}

	for (unsigned int i = 0; i < SPHERE_COUNT; i++)
		world.shapes[i] = makeSphere(&world);
	for (unsigned int i = SPHERE_COUNT; i < SPHERE_COUNT + CUBOID_COUNT; i++)
		world.shapes[i] = makeCuboid(&world);

	/* Plane at the origin */
	world.shapes[SPHERE_COUNT + CUBOID_COUNT - 1] = makeCuboidPlane(&world);

	world.backColor.x = 0.0f;
	world.backColor.y = 0.0f;
	world.backColor.z = 0.0f;
	for (unsigned int i = 0; i < POINT_LIGHT_COUNT; i++)
	{
		world.lights[i] = makePointLight(&world);
		world.backColor.x += world.lights[i].color.x;
		world.backColor.y += world.lights[i].color.y;
		world.backColor.z += world.lights[i].color.z;
	}
	for (unsigned int i = POINT_LIGHT_COUNT;
		i < POINT_LIGHT_COUNT + DISTANT_LIGHT_COUNT; i++)
	{
		world.lights[i] = makeDistantLight(&world);
		world.backColor.x += world.lights[i].color.x;
		world.backColor.y += world.lights[i].color.y;
		world.backColor.z += world.lights[i].color.z;
	}
	for (unsigned int i = POINT_LIGHT_COUNT + DISTANT_LIGHT_COUNT;
		i < POINT_LIGHT_COUNT + DISTANT_LIGHT_COUNT + AREA_LIGHT_COUNT; i++)
	{
		world.lights[i] = makeAreaLight(&world);
		world.backColor.x += world.lights[i].color.x;
		world.backColor.y += world.lights[i].color.y;
		world.backColor.z += world.lights[i].color.z;
	}

	world.backColor.x /= world.lightCount;
	world.backColor.y /= world.lightCount;
	world.backColor.z /= world.lightCount;

	world.fogColor = world.backColor;

	return world;
}

// ---------------------------------
// File to string
// ---------------------------------

char *fileToString(const char *filename)
{
	FILE *f = fopen(filename, "rb");
	if (f)
	{
		fseek(f, 0, SEEK_END);
		unsigned int len = ftell(f);
		fseek(f, 0, SEEK_SET);
		char *buffer = malloc(len + 1);
		if (buffer)
		{
			fread(buffer, 1, len, f);
		}
		fclose(f);
		buffer[len] = '\0';
		return buffer;
	}
	else return NULL;
}

// ---------------------------------
// Warn()
// ---------------------------------

void warn(const char *message) { printf("Warning: %s\n", message); }

// ---------------------------------
// GetBuildReport(), for compilation failure.
// ---------------------------------

const char *getBuildReport(const OpenCL_RT *this, cl_int error)
{
	size_t log_size;
	clGetProgramBuildInfo(this->program, this->device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
	char *log = malloc(log_size);
	clGetProgramBuildInfo(this->program, this->device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
	printf("%s\n", log);
	free(log);
	return "CL_BUILD_PROGRAM_FAILURE\n";
}

/* A nice error message function. */
const char *getErrorString(const OpenCL_RT *this, cl_int error)
{
	switch (error) {
		// run-time and JIT compiler errors
	case 0: return "CL_SUCCESS";
	case -1: return "CL_DEVICE_NOT_FOUND";
	case -2: return "CL_DEVICE_NOT_AVAILABLE";
	case -3: return "CL_COMPILER_NOT_AVAILABLE";
	case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
	case -5: return "CL_OUT_OF_RESOURCES";
	case -6: return "CL_OUT_OF_HOST_MEMORY";
	case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
	case -8: return "CL_MEM_COPY_OVERLAP";
	case -9: return "CL_IMAGE_FORMAT_MISMATCH";
	case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
	case -11: return getBuildReport(this, error);
	case -12: return "CL_MAP_FAILURE";
	case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
	case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
	case -15: return "CL_COMPILE_PROGRAM_FAILURE";
	case -16: return "CL_LINKER_NOT_AVAILABLE";
	case -17: return "CL_LINK_PROGRAM_FAILURE";
	case -18: return "CL_DEVICE_PARTITION_FAILED";
	case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

		// compile-time errors
	case -30: return "CL_INVALID_VALUE";
	case -31: return "CL_INVALID_DEVICE_TYPE";
	case -32: return "CL_INVALID_PLATFORM";
	case -33: return "CL_INVALID_DEVICE";
	case -34: return "CL_INVALID_CONTEXT";
	case -35: return "CL_INVALID_QUEUE_PROPERTIES";
	case -36: return "CL_INVALID_COMMAND_QUEUE";
	case -37: return "CL_INVALID_HOST_PTR";
	case -38: return "CL_INVALID_MEM_OBJECT";
	case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
	case -40: return "CL_INVALID_IMAGE_SIZE";
	case -41: return "CL_INVALID_SAMPLER";
	case -42: return "CL_INVALID_BINARY";
	case -43: return "CL_INVALID_BUILD_OPTIONS";
	case -44: return "CL_INVALID_PROGRAM";
	case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
	case -46: return "CL_INVALID_KERNEL_NAME";
	case -47: return "CL_INVALID_KERNEL_DEFINITION";
	case -48: return "CL_INVALID_KERNEL";
	case -49: return "CL_INVALID_ARG_INDEX";
	case -50: return "CL_INVALID_ARG_VALUE";
	case -51: return "CL_INVALID_ARG_SIZE";
	case -52: return "CL_INVALID_KERNEL_ARGS";
	case -53: return "CL_INVALID_WORK_DIMENSION";
	case -54: return "CL_INVALID_WORK_GROUP_SIZE";
	case -55: return "CL_INVALID_WORK_ITEM_SIZE";
	case -56: return "CL_INVALID_GLOBAL_OFFSET";
	case -57: return "CL_INVALID_EVENT_WAIT_LIST";
	case -58: return "CL_INVALID_EVENT";
	case -59: return "CL_INVALID_OPERATION";
	case -60: return "CL_INVALID_GL_OBJECT";
	case -61: return "CL_INVALID_BUFFER_SIZE";
	case -62: return "CL_INVALID_MIP_LEVEL";
	case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
	case -64: return "CL_INVALID_PROPERTY";
	case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
	case -66: return "CL_INVALID_COMPILER_OPTIONS";
	case -67: return "CL_INVALID_LINKER_OPTIONS";
	case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

		// extension errors
	case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
	case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
	case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
	case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
	case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
	case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
	default: return "Unknown OpenCL error";
	}
}

// ---------------------------------
// Crash()
// ---------------------------------

void crash(const OpenCL_RT *this, const char *message, cl_int status)
{
	printf("\nError: %s, %s", message, getErrorString(this, status));
	system("pause");
	exit(EXIT_FAILURE);
}

// ---------------------------------
// SetPixel()
// -> Faster than SDL_FillRect().
// ---------------------------------

void setPixel(SDL_Surface *dst, unsigned int x, unsigned int y, unsigned int pixel)
{
	int bpp = dst->format->BytesPerPixel;
	unsigned char *p = (unsigned char*)dst->pixels + y * dst->pitch + x * bpp;
	*(unsigned int*)p = pixel;
	return;
	/*
	switch (bpp) {
	case 1:
		*p = (unsigned char)pixel;
		break;
	case 2:
		*(unsigned short*)p = (unsigned short)pixel;
		break;
	case 3:
		if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
			p[0] = (pixel >> 16) & 0xFF;
			p[1] = (pixel >> 8) & 0xFF;
			p[2] = pixel & 0xFF;
		}
		else {
			p[0] = pixel & 0xFF;
			p[1] = (pixel >> 8) & 0xFF;
			p[2] = (pixel >> 16) & 0xFF;
		}
		break;
	case 4:
		*(unsigned int*)p = pixel;
		break;
	}
	*/
}

// ---------------------------------
// RunCLKernel() and Render()
// -> Runs the CL kernel and returns a finished frame.
// ---------------------------------

void runCLKernel(OpenCL_RT *this, SDL_Surface *dst, const World *world)
{
	this->worldBuffer = clCreateBuffer(this->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		sizeof(World), (void*)world, NULL);
	this->toDoFrame = clCreateBuffer(this->context, CL_MEM_WRITE_ONLY,
		(RENDER_WIDTH)*(RENDER_HEIGHT)*sizeof(cl_uint), NULL, NULL);
	this->status = clSetKernelArg(this->kernel, 0, sizeof(cl_mem), (void*)&this->worldBuffer);
	if (this->status != CL_SUCCESS) crash(this, "clSetKernelArg0", this->status);
	this->status = clSetKernelArg(this->kernel, 1, sizeof(cl_mem), (void*)&this->toDoFrame);
	if (this->status != CL_SUCCESS) crash(this, "clSetKernelArg1", this->status);

	this->status = clEnqueueNDRangeKernel(this->commandQueue, this->kernel, 2, NULL,
		this->globalWorkSize, NULL, 0, NULL, NULL);
	if (this->status != CL_SUCCESS) warn("clEnqueueNDRangeKernel");

	this->status = clEnqueueReadBuffer(this->commandQueue, this->toDoFrame, CL_TRUE, 0,
		(RENDER_WIDTH)*(RENDER_HEIGHT)*sizeof(cl_uint), this->doneFrame, 0, NULL, NULL);
	if (this->status != CL_SUCCESS) warn("clEnqueueReadBuffer");

	if (clReleaseMemObject(this->toDoFrame) < 0) warn("clReleaseMemObject0");
	if (clReleaseMemObject(this->worldBuffer) < 0) warn("clReleaseMemObject1");

	int i, j, k, l;
#pragma omp parallel for private(j, k, l)
	for (i = 0; i < RENDER_WIDTH; i++)
	{
		for (j = 0; j < RENDER_HEIGHT; j++)
		{
			for (k = 0; k < MAX_PIXEL_SIZE; k++)
			{
				for (l = 0; l < MAX_PIXEL_SIZE; l++)
				{
					setPixel(dst, i*MAX_PIXEL_SIZE + k, j*MAX_PIXEL_SIZE + l,
						*(this->doneFrame + i + j * RENDER_WIDTH));
				}
			}
		}
	}
}

void render(OpenCL_RT *this, SDL_Surface *dst, const World *world)
{
	runCLKernel(this, dst, world);
	doneRendering = TRUE;
}

// ---------------------------------
// OpenCL initialization
// ---------------------------------

OpenCL_RT initCLKernel()
{
	OpenCL_RT o;
	o.kernelFilename = "ocl_kernel.cl";
	o.kernelString = fileToString(o.kernelFilename);
	o.status = clGetPlatformIDs(1, &o.platform, NULL);
	if (o.status != CL_SUCCESS) crash(&o, "clGetPlatformIDs", o.status);

	o.status = clGetDeviceIDs(o.platform, CL_DEVICE_TYPE_GPU, 1, &o.device, NULL);
	if (o.status != CL_SUCCESS) crash(&o, "clGetDeviceIDs", o.status);

	o.context = clCreateContext(NULL, 1, &o.device, NULL, NULL, NULL);
	o.commandQueue = clCreateCommandQueue(o.context, o.device, 0, NULL);
	o.program = clCreateProgramWithSource(o.context, 1, &o.kernelString, NULL, NULL);

	o.status = clBuildProgram(o.program, 1, &o.device, "-cl-fast-relaxed-math", NULL, NULL);
	if (o.status != CL_SUCCESS) crash(&o, "clBuildProgram", o.status);

	o.globalWorkSize[0] = RENDER_WIDTH;
	o.globalWorkSize[1] = RENDER_HEIGHT;

	o.kernel = clCreateKernel(o.program, "render", NULL);
	return o;
}

void closeCLKernel(OpenCL_RT *o)
{
	if (clReleaseKernel(o->kernel) < 0) warn("clReleaseKernel");
	if (clReleaseProgram(o->program) < 0) warn("clReleaseProgram");
	if (clReleaseCommandQueue(o->commandQueue) < 0) warn("clReleaseCommandQueue");
	if (clReleaseContext(o->context) < 0) warn("clReleaseContext");
	free(o->kernelString);
}

#pragma warning(disable:4100)
int main(int argc, char *argv[])
{
	SDL_Init(SDL_INIT_VIDEO);
	SDL_putenv("SDL_VIDEO_CENTERED=center");
	SDL_WM_SetCaption(SCREEN_TITLE, NULL);
	SDL_Surface *screen = SDL_SetVideoMode(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, SCREEN_FLAGS);
	srand((unsigned int)time(NULL));

	printf("OpenCL C Ray Tracing! Now with stack polymorphism!\n\n");
	printf("Screen resolution: [%d, %d]\n", SCREEN_WIDTH, SCREEN_HEIGHT);
	printf("Render resolution: [%d, %d]\n", RENDER_WIDTH * SUPER_SAMPLES, RENDER_HEIGHT * SUPER_SAMPLES);
	printf("Pixelsize: %d\n", MAX_PIXEL_SIZE);
	printf("Super samples: %dx\n", SUPER_SAMPLES);
	printf("Area light samples: %d\n", LIGHT_SAMPLES);
	printf("Reflection depth: %d\n", MAX_REC_DEPTH);
	printf("Shape count: %d sphere(s), %d cuboid(s)\n",
		SPHERE_COUNT, CUBOID_COUNT);
	printf("Light count: %d point light(s), %d distant light(s), %d area light(s)\n",
		POINT_LIGHT_COUNT, DISTANT_LIGHT_COUNT, AREA_LIGHT_COUNT);

	printf("\n");

	SDL_Event sdl_event;
	World world = makeWorld();
	OpenCL_RT oclrt = initCLKernel();

	unsigned char running = TRUE;
	doneRendering = FALSE;
	while (running)
	{
		size_t startTime = SDL_GetTicks();
		while (SDL_PollEvent(&sdl_event))
		{
			if ((sdl_event.type == SDL_QUIT) ||
				(sdl_event.type == SDL_KEYDOWN &&
				sdl_event.key.keysym.sym == SDLK_ESCAPE)) running = FALSE;
		}

		handleInput(&world);

		if (!doneRendering)
		{
			render(&oclrt, screen, &world);
		}

		SDL_Flip(screen);
		if ((SDL_GetTicks() - startTime) < SCREEN_DELAY)
			SDL_Delay((unsigned int)(SCREEN_DELAY - (SDL_GetTicks() - startTime)));
	}

	closeCLKernel(&oclrt);
	SDL_Quit();
	exit(EXIT_SUCCESS);
}