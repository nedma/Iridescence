/*Iridescence shader made by @xerxes1138, based on the code from : https://belcour.github.io/blog/research/2017/05/01/brdf-thin-film.html
//
//
//
//"A Practical Extension to Microfacet Theory for the Modeling of Varying Iridescence
//Laurent Belcour, Pascal Barla
//ACM Transactions on Graphics (proc. of SIGGRAPH 2017)
//
//May 2017"
//
*/
Shader "Xerxes1138/Iridescence"
{
	Properties
	{
		_Dinc("Dinc", Range(0.0, 10.0)) = 0.570
		_eta2("eta2", Range(1.0, 5.0)) = 1.8
		_eta3("eta3", Range(1.0, 5.0)) = 1.08
		_kappa3("kappa3", Range(0.0, 5.0)) = 0.51
		_roughness("roughness", Range(0.01, 1.0)) = 0.07
		//_IBLTex ("IBL", Cube) = "black" {} // Used to test IS reference, 64 to 128 samples are enough with filtered importance sampling, cube face size is hardcoded at 128px
	}
	SubShader
	{
		Tags { "RenderType"="Opaque" }
		LOD 100

		Pass // One pass forward shader
		{
			Tags {"LightMode"="ForwardBase"}

			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert
			#pragma fragment frag

			#pragma multi_compile_fwdbase nolightmap nodirlightmap nodynlightmap novertexlight
			#pragma multi_compile_fog
			
			#include "UnityCG.cginc"
			#include "Lighting.cginc"
			#include "AutoLight.cginc"

			// Shader parameters
			samplerCUBE _IBLTex;

			half4		_IBLTex_HDR;

			half		_Dinc,	// D is called the Optical Path Difference
						_eta2,	// eta means ratio of indices of refraction, but used as IOR
						_eta3,
						_kappa3,	// k?
						_roughness;

			#include "Iridescence.cginc"


			float3 CalcRi(float2 r12, float2 r23, float2 t121, float OPD, float2 phi2)
			{
				float3 I = 0.0;
				float2 R123 = r12 * r23;
				float2 r123 = sqrt(R123);
				float2 Rs = sqr(t121) * r23 / (1 - R123);

				// Reflectance term for m=0 (DC term amplitude)
				float2 C0 = r12 + Rs;
				float3 S0 = evalSensitivity(0.0, 0.0);
				I += depol(C0) * S0;

				// Reflectance term for m>0 (pairs of diracs)
				float2 Cm = Rs - t121;
				for (int m = 1; m <= 3; ++m)
				{
					Cm *= r123;
					float3 SmS = 2.0 * evalSensitivity(m * OPD, m * phi2.x);
					float3 SmP = 2.0 * evalSensitivity(m * OPD, m * phi2.y);
					I += depolColor(Cm.x * SmS, Cm.y * SmP);
				}

				// Convert back to RGB reflectance
				I = clamp(mul(I, XYZ_TO_RGB), 0.0, 1.0);

				return I;
			}


			// Main function expected by BRDF Explorer
			float3 BRDF(float3 L, float3 V, float3 H, float3 N)
			{
				float Dinc = _Dinc;
				float eta2 = max(_eta2, 1.000277);
				float eta3 = max(_eta3, 1.000277);
				float kappa3 = max(_kappa3, 1e-3);
				float roughness = max(_roughness, 0.05);

				// Force eta_2 -> 1.0 when Dinc -> 0.0
				float eta_2 = lerp(1.0, eta2, smoothstep(0.0, 0.03, Dinc));

				// Compute dot products
				float NdotL = /*saturate*/(dot(N, L));
				float NdotV = /*saturate*/(dot(N, V));
				if (NdotL < 0 || NdotV < 0) 
					return 0.0;

				//H = normalize(L + V);
				float NdotH = /*saturate*/(dot(N, H));
				float VdotH = /*saturate*/(dot(V, H));
				float cosTheta1 = /*saturate*/(dot(H, L));
				float cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1.0 - sqr(cosTheta1)));

				// First interface
				float2 R12, phi12;
				fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
				float2 R21 = R12;
				float2 T121 = 1.0 - R12;
				float2 phi21 = PI - phi12;

				// Second interface
				float2 R23, phi23;
				fresnelConductor(cosTheta2, eta_2, eta3, kappa3, R23, phi23);

				// Phase shift
				float OPD = Dinc * cosTheta2;
				float2 phi2 = phi21 + phi23;

				// Compound terms
				float3 I = CalcRi(R12, R23, T121, OPD, phi2);

				// Microfacet BRDF formula
				float D = GGX(NdotH, roughness);
				float G = smithG_GGX(NdotL, NdotV, roughness);
				return (D*G*I) / (4.0 * NdotL * NdotV);
			}

			float3 Integrate_GGXIridescence(float Roughness, float3 N, float3 V, float3 R, uint NumSamples, int cubeSize )
			{
				float3 SpecularLighting = 0.0;
				#if 1 // No IS
					Unity_GlossyEnvironmentData IBLData;
					IBLData.roughness = Roughness;
					IBLData.reflUVW = R;

					float3 SampleColor = Unity_GlossyEnvironment (UNITY_PASS_TEXCUBE(unity_SpecCube0), unity_SpecCube0_HDR, IBLData);

					float Dinc = _Dinc;
					float eta2 = max(_eta2, 1.000277);
					float eta3 = max(_eta3, 1.000277);
					float kappa3 = max(_kappa3, 1e-3);
					float roughness = max(_roughness, 0.05);

					// Force eta_2 -> 1.0 when Dinc -> 0.0
					float eta_2 = lerp(1.0, eta2, smoothstep(0.0, 0.03, Dinc));

					float cosTheta1 = dot(N, V);

					float cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1.0 - sqr(cosTheta1)));

					// First interface
					float2 R12, phi12;
					fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
					float2 R21 = R12;
					float2 T121 = 1.0 - R12;
					float2 phi21 = PI - phi12;

					// Second interface
					float2 R23, phi23;
					fresnelConductor(cosTheta2, eta_2, eta3, kappa3, R23, phi23);

					// Phase shift
					float OPD = Dinc * cosTheta2;
					float2 phi2 = phi21 + phi23;

					// Compound terms
					float3 I = CalcRi(R12, R23, T121, OPD, phi2);

					SpecularLighting = SampleColor * I;
				#else // IS
					for( uint i = 0; i < NumSamples; i++ )
					{
						float2 Xi = Hammersley2d( i, NumSamples );
			
						float4 H = TangentToWorld(N, ImportanceSampleGGX(Xi, Roughness));

						float3 L = 2.0 * dot( V, H ) * H - V;

						float NoV = saturate( dot( N, V ) );      
	 					float NoL = saturate( dot( N, L ) );
						float NoH = saturate( dot( N, H ) );
						float VoH = saturate( dot( V, H ) );

						if( NoL > 0 )
						{            
							float3 SampleColor = DecodeHDR(texCUBElod(_IBLTex, float4(L, calcLOD(cubeSize, H.w, NumSamples))), _IBLTex_HDR);				
	
							float Dinc = _Dinc;
							float eta2 = max(_eta2, 1.000277);
							float eta3 = max(_eta3, 1.000277);
							float kappa3 = max(_kappa3, 1e-3);
							float roughness = max(_alpha, 0.05);

							// Force eta_2 -> 1.0 when Dinc -> 0.0
							float eta_2 = lerp(1.0, eta2, smoothstep(0.0, 0.03, Dinc));

							// Compute dot products
							float NdotL = dot(N, L);
							float NdotV = dot(N, V);

							//if (NdotL < 0 || NdotV < 0) return 0.0;

							float NdotH = dot(N, H);
							float VdotH = dot(V, H);
	
							float cosTheta1 = dot(H, L);		

							float cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1.0 - sqr(cosTheta1)));

							// First interface
							float2 R12, phi12;
							fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
							float2 R21 = R12;
							float2 T121 = 1.0 - R12;
							float2 phi21 = PI - phi12;

							// Second interface
							float2 R23, phi23;
							fresnelConductor(cosTheta2, eta_2, eta3, kappa3, R23, phi23);

							// Phase shift
							float OPD = Dinc * cosTheta2;
							float2 phi2 = phi21 + phi23;

							// Compound terms
							float3 I = CalcRi(R12, R23, T121, OPD, phi2);

							float G = G_GGX(roughness, NdotL, NdotV);

							SpecularLighting += SampleColor * I * 4.0 * G * NdotL * VdotH / NdotH;	
						}  
					}
					SpecularLighting = SpecularLighting / NumSamples;
				#endif

				return SpecularLighting;
			}

			struct VertexInput
			{
				float4 vertex : POSITION;
				float3 normal : NORMAL;
			};

			struct VertexOutput
			{
				float4 vertex : SV_POSITION;
				float3 worldPos : TEXCOORD0;
				float3 normal : TEXCOORD1;
			};

			VertexOutput vert (VertexInput v)
			{
				VertexOutput o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.worldPos = mul(unity_ObjectToWorld, v.vertex);
				o.normal = UnityObjectToWorldNormal(v.normal.xyz);
				return o;
			}

			float4 frag (VertexOutput i) : SV_Target
			{
				float3 L = _WorldSpaceLightPos0.xyz;
				float3 V = normalize(i.worldPos - _WorldSpaceCameraPos.xyz);
				float3 H = normalize(L + (-V));
				float3 N = normalize(i.normal);
				float3 R = reflect(V, N);

                float3 IBL = Integrate_GGXIridescence(_roughness, N, -V, R, 64 /*NumSamples*/, 128 /*cubeface size*/);
				//IBL = 0;
				float3 brdf = IBL + BRDF(L, -V, H, N) * saturate(dot(N, L)) * _LightColor0.rgb;

				return float4(brdf, 1.0);
			}
			ENDCG
		}
	}
}
