/* Source file for the R3 light class */



/* Include files */

#include "R3Graphics.h"



/* Public variables */

R3PointLight R3null_point_light;



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3PointLight);



/* Public functions */

int 
R3InitPointLight()
{
    /* Return success */
    return TRUE;
}



void 
R3StopPointLight()
{
}



R3PointLight::
R3PointLight(void)
{
}



R3PointLight::
R3PointLight(const R3PointLight& light)
    : R3Light(light),
      position(light.position),
      constant_attenuation(light.constant_attenuation),
      linear_attenuation(light.linear_attenuation),
      quadratic_attenuation(light.quadratic_attenuation)
{
}



R3PointLight::
R3PointLight(const R3Point& position, const RNRgb& color,
	     RNScalar intensity, RNBoolean active,
             RNScalar ca, RNScalar la, RNScalar qa)
    : R3Light(color, intensity, active),
      position(position),
      constant_attenuation(ca),
      linear_attenuation(la),
      quadratic_attenuation(qa)
{
}




void R3PointLight::
SetPosition(const R3Point& position)
{
    // Set position
    this->position = position;
}



void R3PointLight::
SetConstantAttenuation(RNScalar ca)
{
  // Set constant coefficient of attenuation
  this->constant_attenuation = ca;
}



void R3PointLight::
SetLinearAttenuation(RNScalar la)
{
  // Set linear coefficient of attenuation
  this->linear_attenuation = la;
}



void R3PointLight::
SetQuadraticAttenuation(RNScalar qa)
{
  // Set quadratic coefficient of attenuation
  this->quadratic_attenuation = qa;
}



void R3PointLight::
Transform(const R3Transformation& transformation)
{
  // Transform position
  position.Transform(transformation);
}



RNScalar R3PointLight::
IntensityAtPoint(const R3Point& point) const
{
    // Return intensity at point
    if (!IsActive()) return 0.0;
    RNLength d = R3Distance(point, position);
    RNScalar denom = constant_attenuation;
    denom += d * linear_attenuation;
    denom += d * d * quadratic_attenuation;
    if (RNIsZero(denom)) return Intensity();
    else return (Intensity() / denom);
}



R3Vector R3PointLight::
DirectionFromPoint(const R3Point& point) const
{
    // Return direction to point
    R3Vector L = position - point;
    L.Normalize();
    return L;
}



RNScalar R3PointLight::
RadiusOfInfluence(RNScalar intensity_threshhold) const
{
    // Return distance beyond which intensity is below threshold
    // kq*d^2 + kl*d + (kc - 1/a) = 0 (use quadratic formula)
    if (!IsActive()) return 0.0;
    if (RNIsZero(Intensity())) return 0.0;
    if (RNIsZero(intensity_threshhold)) return RN_INFINITY;
    RNScalar A = quadratic_attenuation;
    RNScalar B = linear_attenuation;
    RNScalar C = constant_attenuation - Intensity() / intensity_threshhold;
    RNScalar radius = (-B + sqrt(B*B - 4.0*A*C)) / (2.0*A);
    return radius;
}



R3Sphere R3PointLight::
SphereOfInfluence(RNScalar intensity_threshold) const
{
    // Return sphere within which light intensity is above threshhold
    return R3Sphere(Position(), RadiusOfInfluence(intensity_threshold));
}



RNRgb R3PointLight::
DiffuseReflection(const R3Brdf& brdf, 
    const R3Point& point, const R3Vector& normal) const
{
    // Check if light is active
    if (!IsActive()) return RNblack_rgb;

    // Get material properties
    const RNRgb& Dc = brdf.Diffuse();

    // Get light properties
    const RNRgb& Ic = Color();
    RNScalar I = IntensityAtPoint(point);
    R3Vector L = DirectionFromPoint(point);

    // Compute geometric stuff
    RNScalar NL = normal.Dot(L);
    if (RNIsNegativeOrZero(NL)) return RNblack_rgb;

    // Return diffuse component of reflection
    return (I * NL) * Dc * Ic;
}



RNRgb R3PointLight::
SpecularReflection(const R3Brdf& brdf, const R3Point& eye, 
    const R3Point& point, const R3Vector& normal) const
{
    // Check if light is active
    if (!IsActive()) return RNblack_rgb;

    // Get material properties
    const RNRgb& Sc = brdf.Specular();
    RNScalar s = brdf.Shininess();

    // Get light properties
    const RNRgb& Ic = Color();
    RNScalar I = IntensityAtPoint(point);
    R3Vector L = DirectionFromPoint(point);

    // Compute geometric stuff
    RNScalar NL = normal.Dot(L);
    if (RNIsNegativeOrZero(NL)) return RNblack_rgb;
    R3Vector R = (2.0 * NL) * normal - L;
    R3Vector V = eye - point;
    V.Normalize();
    RNScalar VR = V.Dot(R);
    if (RNIsNegativeOrZero(VR)) return RNblack_rgb;

    // Return specular component of reflection
    return (I * pow(VR,s)) * Sc * Ic;
}



RNRgb R3PointLight::
Reflection(const R3Brdf& brdf, const R3Point& eye, 
    const R3Point& point, const R3Vector& normal) const
{
    // Check if light is active
    if (!IsActive()) return RNblack_rgb;

    // Get material properties
    const RNRgb& Dc = brdf.Diffuse();
    const RNRgb& Sc = brdf.Specular();
    RNScalar s = brdf.Shininess();

    // Get light properties
    RNScalar I = IntensityAtPoint(point);
    R3Vector L = DirectionFromPoint(point);
    const RNRgb& Ic = Color();

    // Compute geometric stuff
    RNScalar NL = normal.Dot(L);
    if (RNIsNegativeOrZero(NL)) return RNblack_rgb;
    R3Vector R = (2.0 * NL) * normal - L;
    R3Vector V = eye - point;
    V.Normalize();
    RNScalar VR = V.Dot(R);

    // Compute diffuse reflection
    RNRgb rgb = (I * NL) * Dc * Ic;

    // Compute specular reflection
    if (RNIsPositive(VR)) rgb += (I * pow(VR,s)) * Sc * Ic;

    // Return total reflection
    return rgb;
}



void R3PointLight::
Draw(int i) const
{
    // Draw light
    GLenum index = (GLenum) (GL_LIGHT0 + i);
    if (index > GL_LIGHT7) return;
    GLfloat buffer[4];
    buffer[0] = Intensity() * Color().R();
    buffer[1] = Intensity() * Color().G();
    buffer[2] = Intensity() * Color().B();
    buffer[3] = 1.0;
    glLightfv(index, GL_DIFFUSE, buffer);
    glLightfv(index, GL_SPECULAR, buffer);
    buffer[0] = Position().X();
    buffer[1] = Position().Y();
    buffer[2] = Position().Z();
    buffer[3] = 1.0;
    glLightfv(index, GL_POSITION, buffer);
    buffer[0] = 180.0;
    glLightf(index, GL_SPOT_CUTOFF, buffer[0]);
    glEnable(index);
    buffer[0] = ConstantAttenuation();
    buffer[1] = LinearAttenuation();
    buffer[2] = QuadraticAttenuation();
    glLightf(index, GL_CONSTANT_ATTENUATION, buffer[0]);
    glLightf(index, GL_LINEAR_ATTENUATION, buffer[1]);
    glLightf(index, GL_QUADRATIC_ATTENUATION, buffer[2]);
}



