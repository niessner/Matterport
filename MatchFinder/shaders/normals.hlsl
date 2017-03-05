
cbuffer ConstantBuffer : register(b0)
{
    matrix worldViewProj; 
    matrix world; 
}


struct VertexShaderOutput
{
	float4 position : SV_POSITION;
	float4 normal : TEXCOORD0;
};

VertexShaderOutput vertexShaderMain(
	float4 position : position,
    float3 normal : normal,
    float4 color : color,
    float2 texCoord : texCoord)
{
	VertexShaderOutput output;
	output.position = mul(position, worldViewProj);
	float3 n = mul(normal, (float3x3)world); //view space normal
	output.normal = float4(n, 0.0f);
	//float4 n = mul(float4(normal, 0.0f), transpose(inverse(world)))
	//output.normal = n;
	return output;
}


float4 pixelShaderMain(VertexShaderOutput input) : SV_Target
{
	//return float4(1.0f, 1.0f, 1.0f, 1.0f);
	return float4(0.5f * (-normalize(input.normal.xyz) + 1.0f), 1.0f); //[-1,1] -> [0,1]
}
