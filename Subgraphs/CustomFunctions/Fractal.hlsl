#ifndef KUYURI_FRACTAL_HLSL
#define KUYURI_FRACTAL_HLSL


//--------------------------------------------------------------------------


#define MAX_RAY_STEPS 64
#define RAY_STOP_TRESHOLD 0.0001
#define MENGER_ITERATIONS 5

float maxcomp(float2 v) { return max(v.x, v.y); }
float3 mod(float3 a, float3 b) { return a - b * floor(a / b); }

float sdCross(float3 p) {
    p = abs(p);
    float3 d = float3(max(p.x, p.y),
                  max(p.y, p.z),
                  max(p.z, p.x));
    return min(d.x, min(d.y, d.z)) - (1.0 / 3.0);
}

float sdCrossRep(float3 p) {
    float3 q = mod(p + 1.0, float3(2.0, 2.0, 2.0)) - 1.0;
    return sdCross(q);
}

float sdCrossRepScale(float3 p, float s) {
    return sdCrossRep(p * s) / s;	
}

float scene(float3 p) {
    float scale = 1.0;
    float dist = 0.0;
    for (int i = 0; i < MENGER_ITERATIONS; i++) {
        dist = max(dist, -sdCrossRepScale(p, scale));
        scale *= 3.0;
    }
    return dist;
}

float3 hsv2rgb(float3 c)
{
    float4 K = float4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    float3 p = abs(frac(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * lerp(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

float4 colorize(float c, float depth) {
	
    float hue = lerp(0.3, 0.9, min(c * 1.2 - 0.05, 1.0));
    float sat = 0.0;
    float lum = lerp(1, 0.0, c);
    float3 hsv = float3(hue, sat, lum);
    float3 rgb = hsv2rgb(hsv);
    return float4(rgb, 1.0);
}


void mainImage_float(float2 fragCoord, float aspectRatio, float3 cameraPos, float3 cameraDir, float fov, out float4 fragColor)
{
    float2 screenPos = fragCoord.xy * 2.0 - 1.0;
	
    float3 cameraPlaneU = float3(1.0, 0.0, 0.0);
    float3 cameraPlaneV = float3(0.0, 1.0, 0.0) * aspectRatio;
    cameraPlaneU *= fov;
    cameraPlaneV *= fov;

    float3 rayPos = cameraPos;
    float3 rayDir = cameraDir + screenPos.x * cameraPlaneU + screenPos.y * cameraPlaneV;
	
    rayDir = normalize(rayDir);
	
    float dist = scene(rayPos);
    int stepsTaken = 0;
    for (int i = 0; i < MAX_RAY_STEPS; i++) {
        if (dist < RAY_STOP_TRESHOLD) {
            continue;
        }
        rayPos += rayDir * dist;
        dist = scene(rayPos);
        stepsTaken = i;
    }
	
    float4 color = colorize(pow(float(stepsTaken) / float(MAX_RAY_STEPS), 1), rayPos.z);
	
    fragColor = color;
}


#endif