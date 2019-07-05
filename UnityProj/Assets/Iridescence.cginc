#pragma once


// Common constants
#define PI 3.14159265358979323846

// XYZ to CIE 1931 RGB color space (using neutral E illuminant)
#define XYZ_TO_RGB float3x3(2.3706743, -0.5138850, 0.0052982,-0.9000405, 1.4253036, -0.0146949, -0.4706338, 0.0885814, 1.0093968)

// Square functions for cleaner code
float sqr(float x) { return x*x; }
float2 sqr(float2 x) { return x*x; }

// Depolarization functions for natural light
float depol(float2 polV) { return 0.5 * (polV.x + polV.y); }
float3 depolColor(float3 colS, float3 colP) { return 0.5 * (colS + colP); }

// GGX distribution function
float GGX(float NdotH, float a)
{
	float a2 = sqr(a);
	float d = sqr(sqr(NdotH) * (a2 - 1.0) + 1.0);
	return a2 / (PI * d + 1e-7);
}

// Smith GGX geometric functions
float smithG1_GGX(float NdotV, float a)
{
	float a2 = sqr(a);
	return 2.0 / (1.0 + sqrt(1.0 + a2 * (1.0 - sqr(NdotV)) / sqr(NdotV)));
}

float smithG_GGX(float NdotL, float NdotV, float a)
{
	return smithG1_GGX(NdotL, a) * smithG1_GGX(NdotV, a);
}

// Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance. Microfacet models forrefraction through rough surfaces. In Proceedings of the 18th Eurographics conference on RenderingTechniques, EGSR'07
float G_GGX(float Roughness, float NdotL, float NdotV)
{
	float m = Roughness;
	float m2 = m * m;

	float G_L = 1.0f / (NdotL + sqrt(m2 + (1 - m2) * NdotL * NdotL));
	float G_V = 1.0f / (NdotV + sqrt(m2 + (1 - m2) * NdotV * NdotV));
	float G = G_L * G_V;

	return G;
}

// Fresnel equations for dielectric/dielectric interfaces.
void fresnelDielectric(in float ct1, in float n1, in float n2, out float2 R, out float2 phi)
{
	float st1 = (1 - ct1 * ct1); // Sinus theta1 'squared'
	float nr = n1 / n2;

	if (sqr(nr) * st1 > 1) // Total reflection
	{
		float2 R = 1.0;

		float2 var = float2(-sqr(nr) *  sqrt(st1 - 1.0 / sqr(nr)) / ct1, -sqrt(st1 - 1.0 / sqr(nr)) / ct1);

		phi = 2.0 * atan(var);
	}
	else // Transmission & Reflection
	{
		float ct2 = sqrt(1 - sqr(nr) * st1);

		float2 r = float2
			(
			(n2 * ct1 - n1 * ct2) / (n2 * ct1 + n1 * ct2),
				(n1 * ct1 - n2 * ct2) / (n1 * ct1 + n2 * ct2)
				);

		phi.x = (r.x < 0.0) ? PI : 0.0;
		phi.y = (r.y < 0.0) ? PI : 0.0;

		R = sqr(r);
	}
}

// Fresnel equations for dielectric/conductor interfaces.
void fresnelConductor(in float ct1, in float n1, in float n2, in float k, out float2 R, out float2 phi)
{
	if (k == 0)// use dielectric formula to avoid numerical issues
	{
		fresnelDielectric(ct1, n1, n2, R, phi);
		return;
	}

	float A = sqr(n2) * (1.0 - sqr(k)) - sqr(n1) * (1.0 - sqr(ct1));
	float B = sqrt(sqr(A) + sqr(2.0 * sqr(n2) * k));
	float U = sqrt((A + B) / 2.0);
	float V = sqrt((B - A) / 2.0);

	R.y = (sqr(n1 * ct1 - U) + sqr(V)) / (sqr(n1 * ct1 + U) + sqr(V));

	float2 var1 = float2(2.0 * n1 * V * ct1, sqr(U) + sqr(V) - sqr(n1* ct1));
	phi.y = atan2(var1.x, var1.y) + PI;

	R.x = (sqr(sqr(n2)*(1 - sqr(k))*ct1 - n1*U) + sqr(2 * sqr(n2)*k*ct1 - n1*V)) / (sqr(sqr(n2)*(1 - sqr(k))*ct1 + n1*U) + sqr(2 * sqr(n2)*k*ct1 + n1*V));

	float2 var2 = float2(2 * n1*sqr(n2)*ct1 * (2 * k*U - (1 - sqr(k))*V), sqr(sqr(n2)*(1 + sqr(k))*ct1) - sqr(n1)*(sqr(U) + sqr(V)));
	phi.x = atan2(var2.x, var2.y);
}

// Evaluation XYZ sensitivity curves in Fourier space
float3 evalSensitivity(float opd, float shift)
{
	// Use Gaussian fits, given by 3 parameters: val, pos and var
	float phase = 2.0 * PI * opd * 1e-6;
	float3 val = float3(5.4856e-13, 4.4201e-13, 5.2481e-13);
	float3 pos = float3(1.6810e+06, 1.7953e+06, 2.2084e+06);
	float3 var = float3(4.3278e+09, 9.3046e+09, 6.6121e+09);
	float3 xyz = val * sqrt(2.0 * PI * var) * cos(pos * phase + shift) * exp(-var * phase * phase);
	xyz.x += 9.7470e-14 * sqrt(2.0 * PI * 4.5282e+09) * cos(2.2399e+06 * phase + shift) * exp(-4.5282e+09 * phase * phase);
	return xyz / 1.0685e-7;
}

float4 TangentToWorld(float3 N, float4 H)
{
	float3 UpVector = abs(N.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
	float3 T = normalize(cross(UpVector, N));
	float3 B = cross(N, T);

	return float4((T * H.x) + (B * H.y) + (N * H.z), H.w);
}

float calcLOD(int cubeSize, float pdf, int NumSamples)
{
	float lod = (0.5 * log2((cubeSize*cubeSize) / float(NumSamples)) + 2.0) - 0.5*log2(pdf);
	return lod;
}

// Brian Karis, Epic Games "Real Shading in Unreal Engine 4"
float4 ImportanceSampleGGX(float2 Xi, float Roughness)
{
	float m = Roughness;
	float m2 = m * m;

	float Phi = 2 * PI * Xi.x;

	float CosTheta = sqrt((1.0 - Xi.y) / (1.0 + (m2 - 1.0) * Xi.y));
	float SinTheta = sqrt(max(1e-5, 1.0 - CosTheta * CosTheta));

	float3 H;
	H.x = SinTheta * cos(Phi);
	H.y = SinTheta * sin(Phi);
	H.z = CosTheta;

	float d = (CosTheta * m2 - CosTheta) * CosTheta + 1;
	float D = m2 / (PI * d * d);
	float pdf = D * CosTheta;

	return float4(H, pdf);
}

uint ReverseBits32(uint bits)
{
#if 0 // Shader model 5
	return reversebits(bits);
#else
	bits = (bits << 16) | (bits >> 16);
	bits = ((bits & 0x00ff00ff) << 8) | ((bits & 0xff00ff00) >> 8);
	bits = ((bits & 0x0f0f0f0f) << 4) | ((bits & 0xf0f0f0f0) >> 4);
	bits = ((bits & 0x33333333) << 2) | ((bits & 0xcccccccc) >> 2);
	bits = ((bits & 0x55555555) << 1) | ((bits & 0xaaaaaaaa) >> 1);
	return bits;
#endif
}

float RadicalInverse_VdC(uint bits)
{
	return float(ReverseBits32(bits)) * 2.3283064365386963e-10; // 0x100000000
}

float2 Hammersley2d(uint i, uint maxSampleCount)
{
	return float2(float(i) / float(maxSampleCount), RadicalInverse_VdC(i));
}







