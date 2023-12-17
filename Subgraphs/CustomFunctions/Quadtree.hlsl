#ifndef KUYURI_TELEPORT_HLSL
#define KUYURI_TELEPORT_HLSL

// https://gist.github.com/patriciogonzalezvivo/670c22f3966e662d2f83
float mod289(float x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
float4 mod289(float4 x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
float4 perm(float4 x){return mod289(((x * 34.0) + 1.0) * x);}

float noise(float3 p){
    float3 a = floor(p);
    float3 d = p - a;
    d = d * d * (3.0 - 2.0 * d);

    float4 b = a.xxyy + float4(0.0, 1.0, 0.0, 1.0);
    float4 k1 = perm(b.xyxy);
    float4 k2 = perm(k1.xyxy + b.zzww);

    float4 c = k2 + a.zzzz;
    float4 k3 = perm(c);
    float4 k4 = perm(c + 1.0);

    float4 o1 = frac(k3 * (1.0 / 41.0));
    float4 o2 = frac(k4 * (1.0 / 41.0));

    float4 o3 = o2 * d.z + o1 * (1.0 - d.z);
    float2 o4 = o3.yw * d.x + o3.xz * (1.0 - d.x);

    return o4.y * d.y + o4.x * (1.0 - d.y);
}

void QuadtreeNoise_float(float2 uv, int iteration, float resolution, float threshold, float3 noiseOffset, float3 noiseScale, out float2 output)
{
    float2 u;
    float2 r = float2(1.0f/resolution, 1.0f/resolution);
    
    float z = r.y;
    for (int i=0; i<iteration; i++) {
        u = floor(uv/z)+.5;
        float2 zu = z*u;
        if (noise(float3(zu.x, zu.y, z) * noiseScale + noiseOffset) < threshold) break;
        z /= 2.;
    }
    output = z/2.-abs(uv-z*u);
}

#endif