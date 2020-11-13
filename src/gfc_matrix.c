#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gfc_matrix.h"
#include "simple_logger.h"

void gfc_matrix_slog(Matrix4 mat)
{
    slog("%f,%f,%f,%f",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
    slog("%f,%f,%f,%f",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
    slog("%f,%f,%f,%f",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
    slog("%f,%f,%f,%f",mat[3][0],mat[3][1],mat[3][2],mat[3][3]);
}

void gfc_matrix_copy(
    Matrix4 d,
    Matrix4 s
  )
{
    if ((!d)||(!s))return;
    if (d == s)return;
    memcpy(d,s,sizeof(Matrix4));
}


void gfc_matrix_multiply(
    Matrix4 out,
    Matrix4 m1,
    Matrix4 m2
  )
{

  out[0][0] = m2[0][0]*m1[0][0] + m2[0][1]*m1[1][0] + m2[0][2]*m1[2][0] + m2[0][3]*m1[3][0];
  out[0][1] = m2[0][0]*m1[0][1] + m2[0][1]*m1[1][1] + m2[0][2]*m1[2][1] + m2[0][3]*m1[3][1];
  out[0][2] = m2[0][0]*m1[0][2] + m2[0][1]*m1[1][2] + m2[0][2]*m1[2][2] + m2[0][3]*m1[3][2];
  out[0][3] = m2[0][0]*m1[0][3] + m2[0][1]*m1[1][3] + m2[0][2]*m1[2][3] + m2[0][3]*m1[3][3];

  out[1][0] = m2[1][0]*m1[0][0] + m2[1][1]*m1[1][0] + m2[1][2]*m1[2][0] + m2[1][3]*m1[3][0];
  out[1][1] = m2[1][0]*m1[0][1] + m2[1][1]*m1[1][1] + m2[1][2]*m1[2][1] + m2[1][3]*m1[3][1];
  out[1][2] = m2[1][0]*m1[0][2] + m2[1][1]*m1[1][2] + m2[1][2]*m1[2][2] + m2[1][3]*m1[3][2];
  out[1][3] = m2[1][0]*m1[0][3] + m2[1][1]*m1[1][3] + m2[1][2]*m1[2][3] + m2[1][3]*m1[3][3];

  out[2][0] = m2[2][0]*m1[0][0] + m2[2][1]*m1[1][0] + m2[2][2]*m1[2][0] + m2[2][3]*m1[3][0];
  out[2][1] = m2[2][0]*m1[0][1] + m2[2][1]*m1[1][1] + m2[2][2]*m1[2][1] + m2[2][3]*m1[3][1];
  out[2][2] = m2[2][0]*m1[0][2] + m2[2][1]*m1[1][2] + m2[2][2]*m1[2][2] + m2[2][3]*m1[3][2];
  out[2][3] = m2[2][0]*m1[0][3] + m2[2][1]*m1[1][3] + m2[2][2]*m1[2][3] + m2[2][3]*m1[3][3];

  out[3][0] = m2[3][0]*m1[0][0] + m2[3][1]*m1[1][0] + m2[3][2]*m1[2][0] + m2[3][3]*m1[3][0];
  out[3][1] = m2[3][0]*m1[0][1] + m2[3][1]*m1[1][1] + m2[3][2]*m1[2][1] + m2[3][3]*m1[3][1];
  out[3][2] = m2[3][0]*m1[0][2] + m2[3][1]*m1[1][2] + m2[3][2]*m1[2][2] + m2[3][3]*m1[3][2];
  out[3][3] = m2[3][0]*m1[0][3] + m2[3][1]*m1[1][3] + m2[3][2]*m1[2][3] + m2[3][3]*m1[3][3];
}

void gfc_matrix_multiply_vector4d(
  Vector4D * out,
  Matrix4    mat,
  Vector4D   vec
)
{
  float x,y,z,w;
  float ox,oy,oz,ow;
  if (!out)return;
  x=vec.x;
  y=vec.y;
  z=vec.z;
  w=vec.w;
  ox=x*mat[0][0] + y*mat[1][0] + mat[2][0]*z + mat[3][0]*w;
  oy=x*mat[0][1] + y*mat[1][1] + mat[2][1]*z + mat[3][1]*w;
  oz=x*mat[0][2] + y*mat[1][2] + mat[2][2]*z + mat[3][2]*w;
  ow=x*mat[0][3] + y*mat[1][3] + mat[2][3]*z + mat[3][3]*w;
  out->x = ox;
  out->y = oy;
  out->z = oz;
  out->w = ow;
}

void gfc_matrix_zero(Matrix4 zero)
{
    memset(zero,0,sizeof(Matrix4));
}

void gfc_matrix_identity(Matrix4 one)
{
    gfc_matrix_zero(one);
    one[0][0] = 1;
    one[1][1] = 1;
    one[2][2] = 1;
    one[3][3] = 1;
}

void gfc_matrix_rotate(
    Matrix4     out,
    Matrix4     m,
    float       degree,
    Vector3D    axis
)
{
    Matrix4 Rotate;
    Matrix4 Result;
    Vector3D temp;
    float a = degree;
    float c = cos(a);
    float s = sin(a);

    //slog("%.2f, %.2f, %.2f", axis.x, axis.y, axis.z);
    vector3d_normalize(&axis);
    
    vector3d_scale(temp,axis,(1 - c));

    Rotate[0][0] = c + temp.x * axis.x;
    Rotate[0][1] = temp.x * axis.y + s * axis.z;
    Rotate[0][2] = temp.x * axis.z - s * axis.y;

    Rotate[1][0] = temp.y * axis.x - s * axis.z;
    Rotate[1][1] = c + temp.y * axis.y;
    Rotate[1][2] = temp.y * axis.z + s * axis.x;

    Rotate[2][0] = temp.z * axis.x + s * axis.y;
    Rotate[2][1] = temp.z * axis.y - s * axis.x;
    Rotate[2][2] = c + temp.z * axis.z;

    Result[0][0] = m[0][0] * Rotate[0][0] + m[1][0] * Rotate[0][1] + m[2][0] * Rotate[0][2];
    Result[0][1] = m[0][1] * Rotate[0][0] + m[1][1] * Rotate[0][1] + m[2][1] * Rotate[0][2];
    Result[0][2] = m[0][2] * Rotate[0][0] + m[1][2] * Rotate[0][1] + m[2][2] * Rotate[0][2];
    Result[0][3] = m[0][3] * Rotate[0][0] + m[1][3] * Rotate[0][1] + m[2][3] * Rotate[0][2];

    Result[1][0] = m[0][0] * Rotate[1][0] + m[1][0] * Rotate[1][1] + m[2][0] * Rotate[1][2];
    Result[1][1] = m[0][1] * Rotate[1][0] + m[1][1] * Rotate[1][1] + m[2][1] * Rotate[1][2];
    Result[1][2] = m[0][2] * Rotate[1][0] + m[1][2] * Rotate[1][1] + m[2][2] * Rotate[1][2];
    Result[1][3] = m[0][3] * Rotate[1][0] + m[1][3] * Rotate[1][1] + m[2][3] * Rotate[1][2];

    Result[2][0] = m[0][0] * Rotate[2][0] + m[1][0] * Rotate[2][1] + m[2][0] * Rotate[2][2];
    Result[2][1] = m[0][1] * Rotate[2][0] + m[1][1] * Rotate[2][1] + m[2][1] * Rotate[2][2];
    Result[2][2] = m[0][2] * Rotate[2][0] + m[1][2] * Rotate[2][1] + m[2][2] * Rotate[2][2];
    Result[2][3] = m[0][3] * Rotate[2][0] + m[1][3] * Rotate[2][1] + m[2][3] * Rotate[2][2];

    Result[3][0] = m[3][0];
    Result[3][1] = m[3][1];
    Result[3][2] = m[3][2];
    Result[3][3] = m[3][3];
    gfc_matrix_copy(out,Result);
}

void gfc_matrix_perspective(
    Matrix4     out,
    float      fov,
    float      aspect,
    float      near,
    float      far
)
{
    float halftanfov = tan(fov * 0.5);
    gfc_matrix_zero(out);

    if (aspect == 0)
    {
        slog("gfc_matrix_perspective: aspect ratio cannot be zero");
        return;
    }
    if (halftanfov == 0)
    {
        slog("gfc_matrix_perspective: bad fov");
        return;
    }
    if (near == far)
    {
        slog("gfc_matrix_perspective: near plane and far plane cannot be the same");
        return;
    }

    gfc_matrix_zero(out);
    out[0][0] = 1 / (aspect * halftanfov);
    out[1][1] = 1 / (halftanfov);
    out[2][2] = - ((far + near) / (far - near));
    out[2][3] = -1;
    if ((far - near) == 0)
    {
        out[3][2] = 0;
    }
    else
    out[3][2] = -(2 * far * near) / (far - near);
    return;
}

void gfc_matrix_view(
    Matrix4  out,
    Vector3D position,
    Vector3D target,
    Vector3D up
)
{
    Vector3D f,s,u;
    vector3d_sub(f,target,position);
    vector3d_normalize(&f);
    
    vector3d_cross_product(&s,f,up);
    vector3d_normalize(&s);
    
    vector3d_cross_product(&u,s,f);
 
    gfc_matrix_identity(out);
    out[0][0] = s.x;
    out[1][0] = s.y;
    out[2][0] = s.z;
    out[0][1] = u.x;
    out[1][1] = u.y;
    out[2][1] = u.z;
    out[0][2] = -f.x;
    out[1][2] = -f.y;
    out[2][2] = -f.z;
    out[3][0] = vector3d_dot_product(s, position)?-vector3d_dot_product(s, position):0;
    out[3][1] = vector3d_dot_product(u, position)?-vector3d_dot_product(u, position):0;
    out[3][2] = vector3d_dot_product(f, position)?vector3d_dot_product(f, position):0;
    
    //gfc_matrix_slog(out);
}

void gfc_matrix_make_translation(
    Matrix4 out,
    Vector3D move
)
{
    if (!out)return;
    gfc_matrix_identity(out);
    out[3][0] = move.x;
    out[3][1] = move.y;
    out[3][2] = move.z;
}

void gfc_matrix_translate(
    Matrix4 out,
    Vector3D move
)
{
    Matrix4 translate,temp;
    gfc_matrix_make_translation(translate,move);
    gfc_matrix_multiply(temp,translate,out);
    gfc_matrix_copy(out,temp);
}

void setRotationX(Matrix4 m_mat, float x)
{
    //gfc_matrix_identity(m_mat);

    m_mat[1][1] = cos(x);
    m_mat[1][2] = sin(x);
    m_mat[2][1] = -sin(x);
    m_mat[2][2] = cos(x);
}

void setRotationY(Matrix4 m_mat, float y)
{
    //gfc_matrix_identity(m_mat);

    m_mat[0][0] = cos(y);
    m_mat[0][2] = -sin(y);
    m_mat[2][0] = sin(y);
    m_mat[2][2] = cos(y);
}

void setRotationZ(Matrix4 m_mat, float z)
{
    //gfc_matrix_identity(m_mat);

    m_mat[0][0] = cos(z);
    m_mat[0][1] = sin(z);
    m_mat[1][0] = -sin(z);
    m_mat[1][1] = cos(z);
}

void setRotation(Matrix4 m_mat, Vector3D rotation) {
    Matrix4 temp;

    Matrix4 rotX;
    gfc_matrix_copy(rotX, m_mat);
    setRotationX(rotX, rotation.x);

    Matrix4 rotY;
    gfc_matrix_identity(rotY);
    setRotationZ(rotY, -rotation.y);

    gfc_matrix_multiply(temp, rotX, rotY);

    Matrix4 rotZ;
    gfc_matrix_identity(rotZ);
    setRotationY(rotZ, rotation.z);

    gfc_matrix_multiply(m_mat, temp, rotZ);
 /*
    double ch = cos(rotation.y);
    double sh = sin(rotation.y);
    double ca = cos(rotation.x);
    double sa = sin(rotation.x);
    double cb = cos(rotation.z);
    double sb = sin(rotation.z);

    m_mat[0][0] = ch * ca;
    m_mat[0][1] = sh * sb - ch * sa * cb;
    m_mat[0][2] = ch * sa * sb + sh * cb;
    m_mat[1][0] = sa;
    m_mat[1][1] = ca * cb;
    m_mat[1][2] = -ca * sb;
    m_mat[2][0] = -sh * ca;
    m_mat[2][1] = sh * sa * cb + ch * sb;
    m_mat[2][2] = -sh * sa * sb + ch * cb;
    */

    /*
    float cosH = cos(rotation.x);
    float sinH = sin(rotation.x);

    float cosP = cos(rotation.y);
    float sinP = sin(rotation.y);

    float cosB = cos(rotation.z);
    float sinB = sin(rotation.z);

    m_mat[0][0] = cosH * cosB + sinH * sinP *sinB;
    m_mat[0][1] = sinB * cosP;
    m_mat[0][2] = -sinH * cosB + cosH * sinP * sinB;
    m_mat[1][0] = -cosH*sinB+sinH*sinP*cosB;
    m_mat[1][1] = cosB * cosP;
    m_mat[1][2] = sinB * sinH + cosH * sinP * cosB;
    m_mat[2][0] = sinH * cosP;
    m_mat[2][1] = -sinP;
    m_mat[2][2] = cosH * cosP;
    */
}   

void setRotation_model(Matrix4 m_mat, Vector3D rotation) {
    Matrix4 temp;

    Matrix4 rotX;
    gfc_matrix_copy(rotX, m_mat);
    setRotationY(rotX, -rotation.y);

    Matrix4 rotY;
    gfc_matrix_identity(rotY);
    setRotationX(rotY, -rotation.x);

    gfc_matrix_multiply(temp, rotX, rotY);

    Matrix4 rotZ;
    gfc_matrix_identity(rotZ);
    setRotationZ(rotZ, rotation.z);

    gfc_matrix_multiply(m_mat, temp, rotZ);
}  
 
float getDeterminant(Matrix4 m_mat)
{
    Vector4D minor, v1, v2, v3;
    float det;

    v1 = vector4d(m_mat[0][0], m_mat[1][0], m_mat[2][0], m_mat[3][0]);
    v2 = vector4d(m_mat[0][1], m_mat[1][1], m_mat[2][1], m_mat[3][1]);
    v3 = vector4d(m_mat[0][2], m_mat[1][2], m_mat[2][2], m_mat[3][2]);

    vector4d_cross_product(&minor, v1, v2, v3);

    det = -(m_mat[0][3] * minor.x + m_mat[1][3] * minor.y + m_mat[2][3] * minor.z +
        m_mat[3][3] * minor.w);
    return det;
}

void gfc_matrix_inverse(Matrix4 out, Matrix4 in)
{
    int a, i, j;
    Matrix4 result;
    Vector4D v, vec[3];
    float det = 0.0f;

    det = getDeterminant(in);
    if (!det) return;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if (j != i)
            {
                a = j;
                if (j > i) a = a - 1;
                vec[a].x = (in[j][0]);
                vec[a].y = (in[j][1]);
                vec[a].z = (in[j][2]);
                vec[a].w = (in[j][3]);
            }
        }
        vector4d_cross_product(&v, vec[0], vec[1], vec[2]);

        //v.cross(vec[0], vec[1], vec[2]);

        result[0][i] = pow(-1.0f, i) * v.x / det;
        result[1][i] = pow(-1.0f, i) * v.y / det;
        result[2][i] = pow(-1.0f, i) * v.z / det;
        result[3][i] = pow(-1.0f, i) * v.w / det;
    }

    gfc_matrix_copy(out, result);
    /*this->setMatrix(out);*/
}


