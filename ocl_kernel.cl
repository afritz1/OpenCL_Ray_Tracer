
/* ocl_kernel.cl; April 2nd, 2015 */

/* Max allowed shapes and lights */
#define MAX_SHAPES 10000
#define MAX_LIGHTS 50

/* Epsilon for vectors */
#define EPSILON 0.001f

/* Size of distant light */
#define DIST_LIGHT_SIZE 0.001f
#define DIST_LIGHT_DOT_MIN (1.0f - DIST_LIGHT_SIZE)

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
	Vector3 point, direction;
	unsigned char depth;
} Ray;

// -----------------------
// Camera typedef
// -----------------------

typedef struct Camera
{
	Vector3 eye, forward, right, up;
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
	Vector3 color;
	unsigned short samples;

	/* PointLight, AreaLight */
	Vector3 point;

	/* DistantLight */
	Vector3 direction;

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
	Vector3 color;

	/* Phong */
	Vector3 specularColor;
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
	Vector3 point;
	Material material;
	
	/* Sphere */
	float radius;
	
	/* Cuboid */
	float width, height, depth;
} Shape;

typedef struct Intersection
{
	float t;
	Vector3 point, normal;
	__global const Shape *shape;
} Intersection;

// -----------------------
// World typedef
// -----------------------

typedef struct World
{
	Camera camera;
	Vector3 backColor, fogColor;	
	unsigned short size, maxRecDepth;
	float fogDist;

	Shape shapes[MAX_SHAPES];
	unsigned short shapeCount;

	Light lights[MAX_LIGHTS];
	unsigned short lightCount;
} World;

// -----------------------
// Camera Methods
// -----------------------

Ray imageRay(__global const Camera *this, float x, float y)
{
	Ray r;
	r.point = this->eye;
	r.direction.x = 
		((this->eye.x + this->forward.x) + (this->up.x - this->right.x) +
		(this->right.x * x * 2.0f) + (-this->up.x * y * 2.0f) - this->eye.x);
	r.direction.y = 
		((this->eye.y + this->forward.y) + (this->up.y - this->right.y) +
		(this->right.y * x * 2.0f) + (-this->up.y * y * 2.0f) - this->eye.y);
	r.direction.z = 
		((this->eye.z + this->forward.z) + (this->up.z - this->right.z) +
		(this->right.z * x * 2.0f) + (-this->up.z * y * 2.0f) - this->eye.z);
	float imgDirLen = 1.0f / native_sqrt(
		(r.direction.x*r.direction.x) +
		(r.direction.y*r.direction.y) +
		(r.direction.z*r.direction.z));
	r.direction.x *= imgDirLen;
	r.direction.y *= imgDirLen;
	r.direction.z *= imgDirLen;
	r.depth = INIT_RAY_DEPTH;
	return r;
}

// -----------------------
// Shape Methods
// -----------------------

Intersection hitSphere(__global const Shape *this, const Ray *ray)
{
	float a =
		(ray->direction.x*ray->direction.x) +
		(ray->direction.y*ray->direction.y) +
		(ray->direction.z*ray->direction.z);
	float b = 2.0f * (
		(ray->direction.x*(ray->point.x - this->point.x)) +
		(ray->direction.y*(ray->point.y - this->point.y)) +
		(ray->direction.z*(ray->point.z - this->point.z)));
	float c =
		(((ray->point.x - this->point.x)*(ray->point.x - this->point.x)) +
		 ((ray->point.y - this->point.y)*(ray->point.y - this->point.y)) +
		 ((ray->point.z - this->point.z)*(ray->point.z - this->point.z))) -
		  (this->radius*this->radius);
	float discriminant = native_sqrt((b * b) - (4.0f * a * c));
	float radiusRecip = 1.0f / this->radius;
	float t1 = (-b + discriminant) * (0.5f * a);
	float t2 = (-b - discriminant) * (0.5f * a);
	Intersection hit;
	hit.t = t2 > 0.0f ? t2 : (t1 > 0.0f ? t1 : FLT_MAX);
	hit.point.x = ray->point.x + (hit.t * ray->direction.x);
	hit.point.y = ray->point.y + (hit.t * ray->direction.y);
	hit.point.z = ray->point.z + (hit.t * ray->direction.z);
	hit.normal.x = (hit.point.x - this->point.x)*radiusRecip;
	hit.normal.y = (hit.point.y - this->point.y)*radiusRecip;
	hit.normal.z = (hit.point.z - this->point.z)*radiusRecip;
	hit.shape = this;
	return hit;
}

Intersection hitCuboid(__global const Shape *this, const Ray *ray)
{
	Intersection intersection;
	intersection.t = MAXFLOAT;
	Vector3 nMin, nMax;
	float tMin, tMax;
	float tX1 = (-this->width + this->point.x - ray->point.x) / ray->direction.x;
	float tX2 = (this->width + this->point.x - ray->point.x) / ray->direction.x;
	if (tX1 < tX2)
	{
		tMin = tX1;
		tMax = tX2;
		nMin.x = -this->width;
		nMin.y = 0.0f;
		nMin.z = 0.0f;
		nMax.x = this->width;
		nMax.y = 0.0f;
		nMax.z = 0.0f;
	}
	else
	{
		tMin = tX2;
		tMax = tX1;
		nMin.x = this->width;
		nMin.y = 0.0f;
		nMin.z = 0.0f;
		nMax.x = -this->width;
		nMax.y = 0.0f;
		nMax.z = 0.0f;
	}

	if (tMin > tMax)
	{
		return intersection;
	}

	float tY1 = (-this->height + this->point.y - ray->point.y) / ray->direction.y;
	float tY2 = (this->height + this->point.y - ray->point.y) / ray->direction.y;

	if (tY1 < tY2)
	{
		if (tY1 > tMin)
		{
			tMin = tY1;
			nMin.x = 0.0f;
			nMin.y = -this->height;
			nMin.z = 0.0f;
		}
		if (tY2 < tMax)
		{
			tMax = tY2;
			nMax.x = 0.0f;
			nMax.y = this->height;
			nMax.z = 0.0f;
		}
	}
	else
	{
		if (tY2 > tMin)
		{
			tMin = tY2;
			nMin.x = 0.0f;
			nMin.y = this->height;
			nMin.z = 0.0f;
		}
		if (tY1 < tMax)
		{
			tMax = tY1;
			nMax.x = 0.0f;
			nMax.y = -this->height;
			nMax.z = 0.0f;
		}
	}

	if (tMin > tMax)
	{
		return intersection;
	}

	float tZ1 = (-this->depth + this->point.z - ray->point.z) / ray->direction.z;
	float tZ2 = (this->depth + this->point.z - ray->point.z) / ray->direction.z;

	if (tZ1 < tZ2)
	{
		if (tZ1 > tMin)
		{
			tMin = tZ1;
			nMin.x = 0.0f;
			nMin.y = 0.0f;
			nMin.z = -this->depth;
		}
		if (tZ2 < tMax)
		{
			tMax = tZ2;
			nMax.x = 0.0f;
			nMax.y = 0.0f;
			nMax.z = this->depth;
		}
	}
	else
	{
		if (tZ2 > tMin)
		{
			tMin = tZ2;
			nMin.x = 0.0f;
			nMin.y = 0.0f;
			nMin.z = this->depth;
		}
		if (tZ1 < tMax)
		{
			tMax = tZ1;
			nMax.x = 0.0f;
			nMax.y = 0.0f;
			nMax.z = -this->depth;
		}
	}

	if (tMin > tMax)
	{
		return intersection;
	}

	if (tMin < 0.0f)
	{
		tMin = tMax;
		nMin = nMax;
	}
	
	if (tMin >= 0.0f)
	{
		intersection.t = tMin;
		intersection.point.x = ray->point.x + (intersection.t * ray->direction.x);
		intersection.point.y = ray->point.y + (intersection.t * ray->direction.y);
		intersection.point.z = ray->point.z + (intersection.t * ray->direction.z);
		float normLen = 1.0f / native_sqrt(
			(nMin.x*nMin.x) + 
			(nMin.y*nMin.y) + 
			(nMin.z*nMin.z));
		intersection.normal.x = nMin.x * normLen;
		intersection.normal.y = nMin.y * normLen;
		intersection.normal.z = nMin.z * normLen;
		intersection.shape = this;
	}

	return intersection;
}

Intersection hitShape(__global const Shape *this, const Ray *ray)
{
	/* This function decides which derived function to call. */
	Intersection noHit;
	switch(this->id)
	{
		case SPHERE_ID:
			return hitSphere(this, ray);
		case CUBOID_ID:
			return hitCuboid(this, ray);
		default:
			noHit.t = MAXFLOAT;
			return noHit;
	}
}

// -----------------------
// Ray Methods
// -----------------------

Vector3 pointAt(const Ray *this, float distance)
{
	Vector3 point;
	point.x = this->point.x + (distance * this->direction.x);
	point.y = this->point.y + (distance * this->direction.y);
	point.z = this->point.z + (distance * this->direction.z);
	return point;
}

Intersection closestHit(const Ray *this, __global const World *world)
{
	Intersection nearestHit;
	nearestHit.t = MAXFLOAT;
	Intersection currentAttempt;

	unsigned int n = world->shapeCount;
	for (unsigned int i = 0; i < n; i++)
	{
		currentAttempt = hitShape(&world->shapes[i], this);

		if (currentAttempt.t < nearestHit.t)
		{
			nearestHit = currentAttempt;
		}
	}

	return nearestHit;
}

// -----------------------
// Light Methods
// -----------------------

Vector3 pointLightDirection(__global const Light *this, const Vector3 *shapePoint)
{
	Vector3 dir;
	dir.x = (this->point.x - shapePoint->x);
	dir.y = (this->point.y - shapePoint->y);
	dir.z = (this->point.z - shapePoint->z);
	float dirLen = 1.0f / native_sqrt(
		(dir.x*dir.x) + 
		(dir.y*dir.y) + 
		(dir.z*dir.z));
	dir.x *= dirLen;
	dir.y *= dirLen;
	dir.z *= dirLen;
	return dir;
}

Vector3 distantLightDirection(__global const Light *this, const Vector3 *shapePoint)
{
	return this->direction;
}

float random()
{
	unsigned int xid = get_global_id(0);
	unsigned int yid = get_global_id(1);
	unsigned int rand = yid * xid * (float)(xid-yid*xid);
	rand *= rand << 32^rand << 16 | rand;
	rand *= rand + (double)(rand);
	return rand;
}

Vector3 areaLightDirection(__global const Light *this, const Vector3 *shapePoint)
{
	Vector3 diff;
	diff.x = (float)(random() * 2) - 1.0f;
	diff.y = (float)(random() * 2) - 1.0f;
	diff.z = (float)(random() * 2) - 1.0f;
	float diffLen = 1.0f / native_sqrt(
		(diff.x*diff.x) + (diff.y*diff.y) + (diff.z*diff.z));
	diff.x *= diffLen * (this->radius * ((float)(random() * 101) / 100.0f));
	diff.y *= diffLen * (this->radius * ((float)(random() * 101) / 100.0f));
	diff.z *= diffLen * (this->radius * ((float)(random() * 101) / 100.0f));
	Vector3 dir =
	{
		(this->point.x + diff.x) - shapePoint->x,
		(this->point.y + diff.y) - shapePoint->y,
		(this->point.z + diff.z) - shapePoint->z
	}; 
	float dirLen = 1.0f / native_sqrt(
		(dir.x*dir.x) + (dir.y*dir.y) + (dir.z*dir.z));
	dir.x *= dirLen;
	dir.y *= dirLen;
	dir.z *= dirLen;
	return dir;
}

Vector3 lightDirection(__global const Light *this, const Vector3 *shapePoint)
{
	/* This function decides which derived function to call. */
	Vector3 noDirection;
	switch(this->id)
	{
		case POINT_LIGHT_ID:
			return pointLightDirection(this, shapePoint);
		case DISTANT_LIGHT_ID:
			return distantLightDirection(this, shapePoint);
		case AREA_LIGHT_ID:
			return areaLightDirection(this, shapePoint);
		default:
			noDirection.x = 0.0f;
			noDirection.y = 0.0f;
			noDirection.z = 0.0f;
			return noDirection;
	}
}

// -----------------------
// Material Methods
// -----------------------

Vector3 flatAt(__global const Material *this, const Intersection *hit, const Ray *ray,
	__global const World *world)
{
	return this->color;
}

Vector3 phongAt(__global const Material *this, const Intersection *hit, const Ray *ray,
	__global const World *world)
{
	/* This method should treat shapes' materials more generically.
	That is, since it must iteratively check the shapes one after
	another, a method more like an inlined materialAt should be used 
	here to allow for reading any shape's material. */



	/* Light polymorphism...? Check ids... */


	
	Material currentMaterial = *this;
	Intersection currentIntersection = *hit;
	Ray currentRay = *ray;

	Vector3 totalColor;
	totalColor.x = 0.0f;
	totalColor.y = 0.0f;
	totalColor.z = 0.0f;

	while (currentIntersection.t < MAXFLOAT && (currentRay.depth <= world->maxRecDepth))
	{
		currentMaterial = currentIntersection.shape->material;
		
		// Cheap placeholder for reflections
		if (currentMaterial.id != PHONG_ID)
		{
			currentIntersection.t = MAXFLOAT;
			break;
		}
		
		Vector3 currentColor;
		currentColor.x = 0.0f;
		currentColor.y = 0.0f;
		currentColor.z = 0.0f;

		Vector3 viewVector;
		viewVector.x = -currentRay.direction.x;
		viewVector.y = -currentRay.direction.y;
		viewVector.z = -currentRay.direction.z;

		float vnDot =
			((viewVector.x * currentIntersection.normal.x) +
			 (viewVector.y * currentIntersection.normal.y) +
			 (viewVector.z * currentIntersection.normal.z));
		char vnSign = vnDot > 0.0f ? 1 : (vnDot < 0.0f ? -1 : 0);

		Vector3 localNormal;
		localNormal.x = currentIntersection.normal.x * vnSign;
		localNormal.y = currentIntersection.normal.y * vnSign;
		localNormal.z = currentIntersection.normal.z * vnSign;

		currentColor.x += currentMaterial.color.x * currentMaterial.ambient;
		currentColor.y += currentMaterial.color.y * currentMaterial.ambient;
		currentColor.z += currentMaterial.color.z * currentMaterial.ambient;

		unsigned int n = world->lightCount;
		for (unsigned int i = 0; i < n; i++)
		{
			__global const Light *light = &world->lights[i];

			Vector3 ambAndDiff;
			ambAndDiff.x = this->color.x * light->color.x;
			ambAndDiff.y = this->color.y * light->color.y;
			ambAndDiff.z = this->color.z * light->color.z;

			Vector3 sampledLight = { 0.0f, 0.0f, 0.0f };
			unsigned int m = light->samples;
			for (unsigned int j = 0; j < m; j++)
			{
				Ray shadowRay;
				shadowRay.direction = lightDirection(light, &currentIntersection.point);
				shadowRay.point.x = currentIntersection.point.x + 
					(shadowRay.direction.x * EPSILON);
				shadowRay.point.y = currentIntersection.point.y + 
					(shadowRay.direction.y * EPSILON);
				shadowRay.point.z = currentIntersection.point.z + 
					(shadowRay.direction.z * EPSILON);

				float lnDot = 
					(shadowRay.direction.x * localNormal.x) + 
					(shadowRay.direction.y * localNormal.y) + 
					(shadowRay.direction.z * localNormal.z);
				char lnSign = lnDot > 0.0f ? 1 : (lnDot < 0.0f ? -1 : 0);

				Vector3 lnReflect;
				lnReflect.x = lnSign * localNormal.x * (2.0f * lnDot) - shadowRay.direction.x;
				lnReflect.y = lnSign * localNormal.y * (2.0f * lnDot) - shadowRay.direction.y;
				lnReflect.z = lnSign * localNormal.z * (2.0f * lnDot) - shadowRay.direction.z;
				float lnReflectLen = 1.0f / native_sqrt(
					(lnReflect.x*lnReflect.x) +
					(lnReflect.y*lnReflect.y) +
					(lnReflect.z*lnReflect.z));
				lnReflect.x *= lnReflectLen;
				lnReflect.y *= lnReflectLen;
				lnReflect.z *= lnReflectLen;

				float lnReflectViewDot = 
					(lnReflect.x * viewVector.x) +
					(lnReflect.y * viewVector.y) +
					(lnReflect.z * viewVector.z);

				Intersection shadowTry = closestHit(&shadowRay, world);
				if (shadowTry.t >= FLT_MAX || !shadowTry.shape->material.shadows)
				{
					// Add up light
					sampledLight.x += ((ambAndDiff.x * this->diffuse) * max(0.0f, lnDot)) +
						((light->color.x * this->specularColor.x) * this->specular *
						native_powr(max(0.0f, lnReflectViewDot), this->shiny));
					sampledLight.y += (ambAndDiff.y * this->diffuse) * max(0.0f, lnDot) +
						((light->color.y  * this->specularColor.y) * this->specular *
						native_powr(max(0.0f, lnReflectViewDot), this->shiny));
					sampledLight.z += (ambAndDiff.z * this->diffuse) * max(0.0f, lnDot) +
						((light->color.z * this->specularColor.z) * this->specular *
						native_powr(max(0.0f, lnReflectViewDot), this->shiny));
				}
			}

			currentColor.x += (sampledLight.x / m);
			currentColor.y += (sampledLight.y / m);
			currentColor.z += (sampledLight.z / m);
		}
		/* Maybe take out these conditionals for now. */
		if (currentRay.depth == INIT_RAY_DEPTH)
		{
			totalColor.x += currentColor.x;
			totalColor.y += currentColor.y;
			totalColor.z += currentColor.z;
		}
		else
		{
			totalColor.x += (currentColor.x * (currentMaterial.specular / currentRay.depth));
			totalColor.y += (currentColor.y * (currentMaterial.specular / currentRay.depth));
			totalColor.z += (currentColor.z * (currentMaterial.specular / currentRay.depth));
		}
		
		/* Attempt another intersection */
		if (currentMaterial.reflective)
		{
			// Make the reflected vector
			Vector3 vnReflect;
			vnReflect.x = vnSign * localNormal.x * (2.0f * vnDot) - viewVector.x;
			vnReflect.y = vnSign * localNormal.y * (2.0f * vnDot) - viewVector.y;
			vnReflect.z = vnSign * localNormal.z * (2.0f * vnDot) - viewVector.z;
			float vnRefLen = 1.0f / native_sqrt(
				(vnReflect.x*vnReflect.x) +
				(vnReflect.y*vnReflect.y) +
				(vnReflect.z*vnReflect.z));
			vnReflect.x *= vnRefLen;
			vnReflect.y *= vnRefLen;
			vnReflect.z *= vnRefLen;

			Ray reflectRay;
			reflectRay.point.x = currentIntersection.point.x + (vnReflect.x * EPSILON);
			reflectRay.point.y = currentIntersection.point.y + (vnReflect.y * EPSILON);
			reflectRay.point.z = currentIntersection.point.z + (vnReflect.z * EPSILON);
			reflectRay.direction = vnReflect;
			reflectRay.depth = currentRay.depth + 1;

			currentIntersection = closestHit(&reflectRay, world);

			// Check if it hits the sky next (because otherwise the background would be missed)
			if (currentIntersection.t == MAXFLOAT)
			{
				totalColor.x += world->fogColor.x * currentMaterial.specular;
				totalColor.y += world->fogColor.y * currentMaterial.specular;
				totalColor.z += world->fogColor.z * currentMaterial.specular;
			}

			currentRay = reflectRay;
		}
		else
		{
			/* This stops the while loop since currentMaterial.reflective is false. */
			currentIntersection.t = MAXFLOAT;
		}
	}

	return totalColor;
}

Vector3 materialAt(__global const Material *this, const Intersection *hit, const Ray *ray,
	__global const World *world)
{
	/* This function decides which derived function to call. */
	Vector3 noColor;
	switch(this->id)
	{
		case FLAT_ID:
			return flatAt(this, hit, ray, world);
		case PHONG_ID:
			return phongAt(this, hit, ray, world);
		default:
			noColor.x = 0.0f;
			noColor.y = 0.0f;
			noColor.z = 0.0f;
			return noColor;
	}	
}

// -----------------------
// World Methods
// -----------------------

Vector3 colorAt(__global const World *this, float x, float y)
{
	Ray ray = imageRay(&this->camera, x, y);
	Intersection nearestHit = closestHit(&ray, this);

	Vector3 color = this->backColor;
	if (nearestHit.t == MAXFLOAT)
	{
		/* Calculate dot between ray and distant light for "sun" effect.
		   In the future, have a "background(Vector3)" function that gets the color
		   of the "skybox", rather than just the lights' colors. */
		for (unsigned int i = 0; i < this->lightCount; i++)
		{
			if (this->lights[i].id != DISTANT_LIGHT_ID) 
				continue;

			Vector3 lightDir = this->lights[i].direction;
			float rayLightDot = 
				(ray.direction.x * lightDir.x) + 
				(ray.direction.y * lightDir.y) + 
				(ray.direction.z * lightDir.z);
			if (rayLightDot > DIST_LIGHT_DOT_MIN)
			{
				Vector3 lightColor = this->lights[i].color;
				color.x += lightColor.x * 1.5f;
				color.y += lightColor.y * 1.5f;
				color.z += lightColor.z * 1.5f;
			}
		}		
	}
	else if (nearestHit.t > this->fogDist)
	{
		color = this->fogColor;
	}
	else
	{
		color = materialAt(
			&nearestHit.shape->material,
			&nearestHit,
			&ray,
			this);
		float fogPortion = nearestHit.t / this->fogDist;
		color.x = (color.x * (1.0f - fogPortion)) + (this->fogColor.x * fogPortion);
		color.y = (color.y * (1.0f - fogPortion)) + (this->fogColor.y * fogPortion);
		color.z = (color.z * (1.0f - fogPortion)) + (this->fogColor.z * fogPortion);
	}

	color.x = color.x > 1.0f ? 1.0f : (color.x < 0.0f ? 0.0f : color.x);
	color.y = color.y > 1.0f ? 1.0f : (color.y < 0.0f ? 0.0f : color.y);
	color.z = color.z > 1.0f ? 1.0f : (color.z < 0.0f ? 0.0f : color.z);

	return color;
}

__kernel void render(__global const World *world, __global unsigned int *output)
{
	unsigned short x = get_global_id(0);
	unsigned short y = get_global_id(1);
	unsigned short width = get_global_size(0);
	unsigned short height = get_global_size(1);
	
	float widthRecip = 1.0f / width;
	float heightRecip = 1.0f / height;
	float superSamplesRecip = 1.0f / world->camera.superSamples;
	float superSamplesSquaredRecip = 1.0f / 
		(world->camera.superSamples * world->camera.superSamples);

	Vector3 color;
	unsigned int colorR = 0;
	unsigned int colorG = 0;
	unsigned int colorB = 0;

	float left = (float)x - 0.5f;
	float right = left + 1.0f;
	float top = (float)y - 0.5f;
	float bottom = top + 1.0f;
	float radius = ((right + bottom - left - top) * 0.25f) * superSamplesRecip;

	for (float i = left + radius; i < right; i += 2.0f * radius)
	{
		for (float j = top + radius; j < bottom; j += 2.0f * radius)
		{
			color = colorAt(world, i * widthRecip, j * heightRecip);
			colorR += (unsigned char)(color.x * 255);
			colorG += (unsigned char)(color.y * 255);
			colorB += (unsigned char)(color.z * 255);
		}
	}

	output[x + y * width] = (unsigned int)
		(((unsigned char)(colorR * superSamplesSquaredRecip) << 16) |
		((unsigned char)(colorG * superSamplesSquaredRecip) << 8) |
		((unsigned char)(colorB * superSamplesSquaredRecip)));
}
