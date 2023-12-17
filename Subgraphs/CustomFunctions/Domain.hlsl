#ifndef KUYURI_DOMAIN_HLSL
#define KUYURI_DOMAIN_HLSL

#define PI 3.1415926535897932384626433832795

/**
 * hash function
 */
uint3 pcg3d( uint3 s ) {
  s = s * 1145141919u + 1919810u;
  s.x += s.y * s.z;
  s.y += s.z * s.x;
  s.z += s.x * s.y;
  s ^= s >> 16;
  s.x += s.y * s.z;
  s.y += s.z * s.x;
  s.z += s.x * s.y;
  return s;
}

/**
 * pcg3d but float
 */
float3 pcg3df( float3 s ) {
  uint3 r = pcg3d( asuint( s ) );
  return float3( r ) / float( 0xffffffffu );
}

struct GridResult {
  float3 cell;
  float d;
  float3 n;
};

/**
 * Return the grid cell information the ray currently belongs to
 * and a distance to the boundary of the cell.
 */
GridResult gridTraversal( float3 ro, float3 rd ) {
  GridResult result;

  result.cell = floor( ro + rd * 1E-3 ) + 0.5;
  
  // calculate the distance to the boundary
  // It's basically a backface only cube intersection
  // See the iq shader: https://www.shadertoy.com/view/ld23DV
  float3 src = -( ro - result.cell ) / rd;
  float3 dst = abs( 0.5 / rd );
  float3 bv = src + dst;
  result.d = min( min( bv.x, bv.y ), bv.z );

  result.n = -step( bv, float3( result.d, result.d, result.d ) ) * sign( rd );
  
  return result;
}

float random(float v)
{
  return frac(sin(691.43 * v) * 571.54);
}

float3x3 rotate3d(float3 xyz)
{
  float a = xyz.y;
  float b = xyz.x;
  float c = xyz.z;
  return float3x3(cos(b)*cos(c), sin(a)*sin(b)*cos(c) - cos(a)*sin(c), cos(a)*sin(b)*cos(c) + sin(a)*sin(c),
              cos(b)*sin(c), sin(a)*sin(b)*sin(c) + cos(a)*cos(c), cos(a)*sin(b)*sin(c) - sin(a)*cos(c),
              -sin(b),       sin(a)*cos(b),                        cos(a)*cos(b));
}

float sdsphere( float3 p, float r ) {
  return length( p ) - r;
}

float sdbox( float3 p, float3 s ) {
  float3 d = abs( p ) - s;
  return length( max( d, 0.0 ) ) + min( 0.0, max( max( d.x, d.y ), d.z ) );
}

float sdtorus( float3 p, float r, float R ) {
  float2 d = float2( length( p.xy ) - R, p.z );
  return length( d ) - r;
}

float sdBoxFrame( float3 p, float3 b, float e )
{
  p = abs(p  )-b;
  float3 q = abs(p+e)-e;
  return min(min(
      length(max(float3(p.x,q.y,q.z),0.0))+min(max(p.x,max(q.y,q.z)),0.0),
      length(max(float3(q.x,p.y,q.z),0.0))+min(max(q.x,max(p.y,q.z)),0.0)),
      length(max(float3(q.x,q.y,p.z),0.0))+min(max(q.x,max(q.y,p.z)),0.0));
}

float mod( float a, float b ) {
  return a - b * floor( a / b );
}

float3 mod( float3 a, float3 b ) {
  return a - b * floor( a / b );
}

float2 path(float z) {
  return float2(0.6, 0.0);
}


float map( float3 p, GridResult grid, float iTime ) {
  float3 gridCell = grid.cell;
  
  float3 po = p;
  p -= gridCell;
  
  float kind = mod( dot( gridCell, float3( 1.0, 1.0, 1.0 ) ), 3.0 );

  float3 p0 = float3(0, 0, 0);
  float3 p1 = float3(
    random(gridCell.x * 215 + gridCell.y * 35 + gridCell.z * 784) * 2.0 - 1.0,
    random(gridCell.x * 436 + gridCell.y * 55 + gridCell.z * 885) * 2.0 - 1.0,
    random(gridCell.x * 375 + gridCell.y * 886 + gridCell.z * 444) * 2.0 - 1.0
  );
  p1 *= 0.2;
  float3 move = lerp(p0, p1, smoothstep(0.4, 0.6, cos(iTime * 4.0) * 0.5 + 0.5));

  float3 posRandom = float3(
    random(gridCell.x * 895 + gridCell.y * 860 + gridCell.z * 34),
    random(gridCell.x * 223 + gridCell.y * 382 + gridCell.z * 967),
    random(gridCell.x * 116 + gridCell.y * 647 + gridCell.z * 38)
  );
  posRandom = posRandom * 2.0 - 1.0;
  posRandom *= 0.4;
  
  float3 sizeRandom = float3(
    random(gridCell.x * 27 + gridCell.y * 49 + gridCell.z * 2),
    random(gridCell.x * 990 + gridCell.y * 85 + gridCell.z * 77),
    random(gridCell.x * 563 + gridCell.y * 267 + gridCell.z * 85)
  );
  //float3 sizeRandom = random(gridCell.x * 27 + gridCell.y * 49 + gridCell.z * 2);
  sizeRandom = sizeRandom.x * sizeRandom.x * sizeRandom.x;
  float sizeMin = 0.1;
  float sizeMax = 1;
  sizeRandom = lerp(float3(sizeMin, sizeMin, sizeMin), float3(sizeMax, sizeMax, sizeMax), sizeRandom);
  sizeRandom = 0.3;

  float d = 1E-2;
  //p += posRandom;
  p += posRandom;
  d = sdbox(p,sizeRandom);

  for(int i = 0; i < 20; i++)
  {
    float3 pr = float3(
      random(gridCell.x * 895 + gridCell.y * 860 + gridCell.z * 34 + i * 1),
      random(gridCell.x * 223 + gridCell.y * 382 + gridCell.z * 967 + i * 4),
      random(gridCell.x * 116 + gridCell.y * 647 + gridCell.z * 38 + i * 50)
    );
    float sizeRatio = i < 10 ? 0.5 : 0.3;
    pr = pr * 2.0 - 1.0;
    pr *= 0.3;
    d = min(d, sdbox(p + pr, sizeRandom * sizeRatio));
  }

  d = max(-sdbox(po, float3(1, 1, 1000)), d);

  return d;
}

float3 nmap( float3 p, GridResult grid, float iTime ) {
  float2 d = float2( 0, 1E-4 );
  return normalize( float3(
    map( p + d.yxx, grid, iTime ) - map( p - d.yxx, grid, iTime ),
    map( p + d.xyx, grid, iTime ) - map( p - d.xyx, grid, iTime ),
    map( p + d.xxy, grid, iTime ) - map( p - d.xxy, grid, iTime )
  ) );
}

float Outline_OffsetDifference(float3 p, GridResult grid, float iTime)
{
  // 0.003 is the offset size, and thus outline thickness, in uv
  float2 d = float2( 0, 1.0 * 1E-2 );         

  float2 marchOA = nmap( p - d.yxx, grid, iTime );
  float2 marchOB = nmap( p - d.xyx, grid, iTime );
  float2 marchA = nmap( p + d.yxx, grid, iTime );
  float2 marchB = nmap( p + d.xyx, grid, iTime );
    
  // 0.07 is the depth threshold is world units, and thus is dependent on scene geometry for a proper value.
  float diff = clamp(max(length(marchOA - marchA), length(marchOB - marchB)) / 0.07, 0.0, 1.0);
    
  // 0.6 is a control value for outline stroke thickness, and 8.0 is stroke strength.
  return 1.0 - smoothstep(0.6, -0.001, pow(diff, 10.0));
}

float shadow(float3 ro,float3 rd, GridResult grid, float iTime){
  float t = 0.;
  for(int i=0;i<200;i++){
    float d = map(ro+rd*t, grid, iTime);
    if(abs(d)<0.001){
      return 0.;
    }
    t+=d;
  }
    
  return 1.;
}

inline float3 GetCameraPosition()    { return UNITY_MATRIX_I_V._m03_m13_m23; }
inline float3 GetCameraForward()     { return -UNITY_MATRIX_V[2].xyz;    }
inline float3 GetCameraUp()          { return UNITY_MATRIX_V[1].xyz;     }
inline float3 GetCameraRight()       { return UNITY_MATRIX_V[0].xyz;     }
inline float  GetCameraFocalLength() { return abs(UNITY_MATRIX_P[1][1]); }
inline float  GetCameraNearClip()    { return _ProjectionParams.y;       }
inline float  GetCameraFarClip()     { return _ProjectionParams.z;       }
inline bool   IsCameraPerspective()  { return any(UNITY_MATRIX_P[3].xyz); }
inline bool   IsCameraOrtho()        { return !IsCameraPerspective(); }

inline float3 _GetCameraDirection(float2 sp)
{
  float3 camDir      = GetCameraForward();
  float3 camUp       = GetCameraUp();
  float3 camSide     = GetCameraRight();
  float  focalLen    = GetCameraFocalLength();
  return normalize((camSide * sp.x) + (camUp * sp.y) + (camDir * focalLen));
}

inline float3 GetCameraDirection(float4 projPos)
{
  projPos.xy /= projPos.w;
  projPos.xy = (projPos.xy - 0.5) * 2.0;
  projPos.x *= _ScreenParams.x / _ScreenParams.y;
  return _GetCameraDirection(projPos.xy);
}

inline float GetDistanceFromCameraToNearClipPlane(float4 projPos)
{
  projPos.xy /= projPos.w;
  projPos.xy = (projPos.xy - 0.5) * 2.0;
  projPos.x *= _ScreenParams.x / _ScreenParams.y;
  float3 norm = normalize(float3(projPos.xy, GetCameraFocalLength()));
  return GetCameraNearClip() / norm.z;
}

void mainImage_float(float2 fragCoord, float aspectRatio, float3 ro, float3 cd, float3 cuu, float fov, float iTime, float3 centerPos, float centerRange, float3 edgeColor, float3 col1, float3 col2, float3 skyCol, out float4 fragColor){
  float2 uv = fragCoord;
  float2 p = uv * 2.0 - 1.0;
  p.x *= aspectRatio;

  float3 cs = normalize(cross(cd, cuu));
  float3 cu = normalize(cross(cs, cd));
  float td = 1.0 / tan(fov / 2.0);
  float3 rd = normalize(cs * p.x + cu * p.y + cd * td);
  
  float rl = 1E-2;
  float3 rp = ro + rd * rl;
  float gridlen = 0.0;
  GridResult grid;
  float dist;
  float rlMin = 1E10;
  float ac = 0.0;
  float3 g = float3( 0.0, 0.0, 0.0 );
  
  for( int i = 0; i < 64; i ++ ) {
    if ( gridlen <= rl ) {
      grid = gridTraversal( rp, rd );
      gridlen += grid.d;
    }

    float d = map( rp, grid, iTime );
    
    dist = d;
    //dist = max(abs(dist), 1e-3);
    rl = min( rl + dist, gridlen );
    rp = ro + rd * rl;

    d = map( rp.xyz, grid, iTime ) + map( rp.yzx, grid, iTime ) + map( rp.zxy, grid, iTime );
    g += float3(d, d*d, d*d*d);

    ac += exp(-dist * 3.0);
    rlMin = min(rlMin, rl);
  }
  ac *= 0.01;

  float3 colAl = lerp(col2, col1, saturate(distance(rp, centerPos) / centerRange));
  float3 col = colAl;
  dist = map( rp, grid, iTime );
  float edge = saturate(Outline_OffsetDifference( rp, grid, iTime));
  col += edge * edgeColor;
  float3 emit = pow(dist + 2.0, -1.2) * colAl;
  col += emit;

  float3 lightdir = normalize(float3(1.,2.,-1.));
  //col += 1.0 - shadow(rp+0.01*n,lightdir, grid.cell, iTime);
  float3 fog = lerp(0, 1, exp( -0.08 * rl )); // distance fog
  ac = ac == 1 ? 1 : 1 - pow(2, -10 * ac);
  col = lerp(skyCol, col * ac, 1 - pow(1 - ac, 3));
  //col *= lerp(0.0, 1.0, exp( -0.1 * rl )); // distance fog
  fragColor = float4( col, 1.0 );
}

#endif