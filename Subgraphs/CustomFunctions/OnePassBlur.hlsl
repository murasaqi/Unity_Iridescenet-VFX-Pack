
void OnePassBlur_float(Texture2D tex, SamplerState texSampler, float2 uv, float2 direction, float radius, out float4 outColor)
{
    float Pi = 3.1415926535897932384626433832795f;
    float dirLength = 30.0f;
    float quality = 50.0f;
    
    float2 rad = radius * direction;

    float4 color = tex.Sample(texSampler, uv);
    for(float d = 0.0f; d < Pi; d += Pi / dirLength)
    {
        for(float i = 1.0f / quality; i <= 1.0f; i += 1.0f / quality)
        {
            float2 offset = float2(cos(d), 0) * i * rad;
            color += tex.Sample(texSampler, saturate(uv + offset));
        }
    }

    color /= quality * dirLength - 15.0f;
    outColor = color;
}