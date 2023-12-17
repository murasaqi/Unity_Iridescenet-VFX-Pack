#ifndef KUYURI_TUNNEL_HLSL
#define KUYURI_TUNNEL_HLSL

float hash(float n) {
	return frac(sin(n)*43758.5453);
}

float3 hash(float3 p){ 
    
    float n = sin(dot(p, float3(7, 157, 113)));    
    return frac(float3(2097152, 262144, 32768)*n); 
}

float noise(float g) {
	float p = floor(g);
	float f = frac(g);

	return lerp(hash(p), hash(p + 1.0), f);
}

float voronoi(float3 x) {
	float3 p = floor(x);
	float3 f = frac(x);

	float2 res = float2(8.0, 8.0);

	for(int i = -1; i <= 1; i++)
	for(int j = -1; j <= 1; j++)
	for(int k = -1; k <= 1; k++) {
		float3 g = float3(float(i), float(j), float(k));
		float3 r = g + hash(p + g) - f;

		float d = max(abs(r.x), max(abs(r.y), abs(r.z)));

		if(d < res.x) {
			res.y = res.x;
			res.x = d;
		} else if(d < res.y) {
			res.y = d;
		}
	}

	return res.y - res.x;
}

float2 path(float z) {
	return float2(12.0, 0.0);
}

float map(float3 p) {
	float4 q = float4(p, 1.0);
	q.x += 1.0;
	float2 tun = abs(p.xy - path(p.z))*float2(0.6, 0.3);

	return min(0.25*abs(q.y)/q.w, 1.0 - max(tun.x, tun.y));
}

float march(float3 ro, float3 rd, float mx) {
	float t = 0.0;

	for(int i = 0; i < 200; i++) {
		float d = map(ro + rd*t);
		if(d < 0.001 || t >= mx) break;
		t += d*0.75;
	}

	return t;
}

float3 normal(float3 p) {
	float2 h = float2(0.001, 0.0);
	float3 n = float3(
		map(p + h.xyy) - map(p - h.xyy),
		map(p + h.yxy) - map(p - h.yxy),
		map(p + h.yyx) - map(p - h.yyx)
	);
	return normalize(n);
}

float ao(float3 p, float3 n) {
    float o = 0.0, s = 0.005;
    for(int i = 0; i< 15; i++) {
        float d = map(p + n*s);
        o += (s - d);
        s += s/float(i + 1);
    }
    
    return 1.0 - clamp(o, 0.0, 1.0);
}

float3 material(float3 p, float iTime) {
	p.z *= 0.01;
	float v = 0.0;
	float a = 0.7, f = 1.0;

	for(int i = 0; i < 4; i++) {
		float v1 = voronoi(p*f + 5.0);
		float v2 = 0.0;

		if(i > 0) {
			v2 = voronoi(p*f*0.1 + 50.0 + 0.15*iTime);

			float va = 0.0, vb = 0.0;
			va = 1.0 - smoothstep(0.0, 0.1, v1);
			vb = 1.0 - smoothstep(0.0, 0.08, v2);
			v += a*pow(va*(0.5 + vb), 4.0);
		}

		v1 = 1.0 - smoothstep(0.0, 0.3, v1);
		v2 =  a*noise(v1*5.5 + 0.1);

		v += v2;

		f *= 3.0;
		a *= 0.5;
	}

	return float3(pow(v, 6.0), pow(v, 4.0), pow(v, 2.0))*2.0;
}

float3x3 camera(float3 eye, float3 lat) {
	float3 ww = normalize(lat - eye);
	float3 uu = normalize(cross(float3(0, 1, 0), ww));
	float3 vv = normalize(cross(ww, uu));

	return float3x3(uu, vv, ww);
}

void mainImage_float(float2 fragCoord, float aspectRatio, float iTime, float3 ro, float3 la, float fov, out float4 fragColor) {
	float2 uv = -1.0 + 2.0*fragCoord.xy;
	uv.x *= aspectRatio;

	float3 col = float3(0, 0, 0);

	ro.xy += path(ro.z);
	la.xy += path(la.z);
	float3 rd = normalize(mul(camera(ro, la), float3(uv, fov)));

	const float tmax = 20.0;
	float i = march(ro, rd, tmax);
	if(i < tmax) {
		float3 pos = ro + rd*i;
		float3 nor = normal(pos);

		float3 rig = ro + float3(0, 0, 3);
		rig.xy += path(rig.z);
		float3 key = normalize(pos - rig);

		col  = 0.1*float3(0, 0, 1);
		col += 0.9*clamp(dot(key, nor), 0.0, 1.0)*float3(1.0/max(1.0, i), 1, 1);
		col += 0.4*clamp(dot(-key, nor), 0.0, 1.0)*float3(1.0/max(1.0, i), 1, 1);

		col *= material(pos, iTime);
	}

	col = lerp(col, float3(0, 0, 0), 1.0 - exp(-0.6*i));
    
    col = 1.0 - exp(-0.5*col);
	float tmp = 1.0/2.2;
    col = pow(abs(col), float3(tmp, tmp, tmp));

	fragColor = float4(col, 1);
}

#endif