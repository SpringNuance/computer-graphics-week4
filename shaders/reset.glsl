#ifdef ComputeShader

#define SOLUTION // tag for solution parser, ignore
layout (local_size_x = 16, local_size_y = 16) in;

struct clothState
{
	vec3 pos, vel;
};

layout (std430, binding=0) buffer cloth_state
{ 
	clothState state[];
};

layout(std140, binding=3) uniform params
{
	float springK;
	float springDamp;
	float dragK;
	float windStrength;
	float springBreakThreshold;
	float scale;
	float dt;

	// Bullet projectile data
	float bulletR;
	vec4 bulletPos;
	vec4 bulletVel;
	
	int w, h;
};

void main() 
{
	int idx = int(gl_GlobalInvocationID.y * w + gl_GlobalInvocationID.x);

	if (gl_GlobalInvocationID.y < h && gl_GlobalInvocationID.x < w)
	{
		// YOUR CODE HERE (compute extra):
		// Write the initial state of the cloth system into the state buffer.

		// state[idx].pos = ...
		// state[idx].vel = ...

#ifdef SOLUTION
		vec2 gridPos = vec2(gl_GlobalInvocationID.xy) / vec2(w, h) * 2.0 - vec2(1.0);

		state[idx].pos = vec3(gridPos.x, 1.0, 1.0f + gridPos.y);
		state[idx].vel = vec3(.0);
#endif
	}
}

#endif