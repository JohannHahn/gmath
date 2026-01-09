#ifndef GMATH_HPP
#define GMATH_HPP

#include <cstring>
#include <string>
#include <print>
#include <assert.h>
#include <math.h>
#include <ostream>

namespace gmath {

struct Mat4;
struct Mat3;

constexpr float PI = 3.14159265359f;

struct Vec3 {
    float x;
    float y;
    float z;
	
    float length() const {
        return sqrtf(x*x + y*y + z*z);
    }

    std::string to_str() const {
	return std::string("Vec3: ") + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z);
    }

    void multiply(const Mat4& m);
    void multiply(const Mat3& m);

    Vec3 project(int width, int height) const {
	Vec3 v;
	v.z = z;
	if (v.z == 0.f) v.z = 0.00001f;

	//std::println("projecting {} to width = {}, height = {}", to_str(), width, height);
	float clip_x = x / v.z;
	float clip_y = y / v.z;
	//std::println("clipspace x = {}, y = {}", clip_x, clip_y);
	
	// move zero to top left and scale 
	v.x = ((clip_x + 1.f) / 2.f) * width;
	v.y = ((-clip_y + 1.f ) / 2.f) * height;

	//std::println("final form = {}", v.to_str());

	return v;
    }
};


struct Vec4 {
    float x;
    float y;
    float z;
    float w;
	
    float length() const {
        return sqrtf(x*x + y*y + z*z + w*w);
    }

    std::string to_str() {
	return std::string("Vec3: ") + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ", " + std::to_string(w);
    }

    void multiply(const Mat4& m);

};

struct Mat2 {
    float data[2*2] = {0.f};
    int size = 2;

    void multiply(float a) {
	for (float& f : data) {
	    f *= a;
	}
    }

    void multiply(const Mat2& m) {
	float new_data[2*2] = {0};

	for (int row = 0; row < size; ++row) {
	    int start = row * size;

	    for (int col = 0; col < size; ++col) {
		new_data[col + row * size] = data[start] * m.data[col];

		for (int i = 1; i < size; ++i) {
		    new_data[col + row * size] += data[start + i] * m.data[col + size * i];
		}
	    }
	    
	}
	memcpy(data, new_data, sizeof(float) * size * size);
    }

    void zero() {
	std::memset(data, 0.f, sizeof(float) * size * size);
    }

    std::string to_str() {
	std::string out = "";
	for (int y = 0; y < size; ++y) {
	    for (int x = 0; x < size; ++x) {
		out += std::to_string(data[x + y * size]) + " ";
	    }
	    out += "\n";
	}
	return out;
    }

    static Mat2 identity() {
	return {{1.f, 0.f,  0.f, 1.f}};
    } 
};

struct Mat3 {
    float data[3*3] = {0.f};
    int size = 3;

    void multiply(float a) {
	for (float& f : data) {
	    f *= a;
	}
    }
    void multiply(const Mat3& m) {
	float new_data[3*3] = {0};

	for (int row = 0; row < size; ++row) {
	    int start = row * size;

	    for (int col = 0; col < size; ++col) {
		new_data[col + row * size] = data[start] * m.data[col];

		for (int i = 1; i < size; ++i) {
		    new_data[col + row * size] += data[start + i] * m.data[col + size * i];
		}
	    }
	    
	}
	memcpy(data, new_data, sizeof(float) * size * size);
    }

    void zero() {
	std::memset(data, 0.f, sizeof(float) * size * size);
    }

    std::string to_str() {
	std::string out = "";
	for (int y = 0; y < size; ++y) {
	    for (int x = 0; x < size; ++x) {
		out += std::to_string(data[x + y * size]) + " ";
	    }
	    out += "\n";
	}
	return out;
    }

    static Mat3 identity() {
	return {{1.f, 0.f, 0.f,  0.f, 1.f, 0.f,  0.f, 0.f, 1.f}};
    } 
};

Mat4 operator*(const Mat4& a, const Mat4& b);

struct Mat4 {
    float data[4*4] = {0.f};
    int size = 4;

    void multiply(float a) {
	for (float& f : data) {
	    f *= a;
	}
    }
    void multiply(const Mat4& m) {
	float new_data[4*4] = {0};

	for (int row = 0; row < size; ++row) {
	    int start = row * size;

	    for (int col = 0; col < size; ++col) {
		new_data[col + row * size] = data[start] * m.data[col];

		for (int i = 1; i < size; ++i) {
		    new_data[col + row * size] += data[start + i] * m.data[col + size * i];
		}
	    }
	    
	}
	memcpy(data, new_data, sizeof(float) * size * size);
    }

    void zero() {
	std::memset(data, 0.f, sizeof(float) * size * size);
    }

    std::string to_str() {
	std::string out = "";
	for (int y = 0; y < size; ++y) {
	    for (int x = 0; x < size; ++x) {
		out += std::to_string(data[x + y * size]) + " ";
	    }
	    out += "\n";
	}
	return out;
    }

    static Mat4 identity() {
	return {{1.f, 0.f, 0.f, 0.f,  
		 0.f, 1.f, 0.f, 0.f,  
		 0.f, 0.f, 1.f, 0.f,  
	         0.f, 0.f, 0.f, 1.f}};
    } 

    static Mat4 translation(const Vec3& v) {
	return {{1.f, 0.f, 0.f, v.x,  
		 0.f, 1.f, 0.f, v.y,  
		 0.f, 0.f, 1.f, v.z,  
		 0.f, 0.f, 0.f, 1.f}};
    }

    static Mat4 rotation_x(float th) {
	return {{1.f, 0.f, 0.f, 0.f,  
		 0.f, cos(th), -sin(th), 0.f,  
		 0.f, sin(th), cos(th), 0.f,  
		 0.f, 0.f, 0.f, 1.f}};
    }

    static Mat4 rotation_y(float th) {
	return {{cos(th), 0.f, sin(th), 0.f,  
		0.f,	  1.f, 0.f,     0.f,  
		-sin(th), 0.f, cos(th), 0.f,  
		 0.f,     0.f, 0.f,     1.f}};
    }

    static Mat4 rotation_z(float th) {
	return {{cos(th), -sin(th), 0.f, 0.f, 
		sin(th), cos(th),   0.f, 0.f,  
		 0.f, 0.f, 1.f, 0.f,  
	         0.f, 0.f, 0.f, 1.f}};
    }

    static Mat4 get_model(Vec3 pos, Vec3 angles) {
	Mat4 m = identity();
	m = translation(pos);
	Mat4 rotation = rotation_x(angles.x) * rotation_y(angles.y) * rotation_z(angles.z);
	m.multiply(rotation);
	return m;
    }

    static Mat4 inverse(const Mat4 m) {

	Mat4 r;

	// Transponiere 3x3 Rotation
	r.data[0] = m.data[0]; r.data[1] = m.data[4]; r.data[2] = r.data[8];
	r.data[4] = m.data[1]; r.data[5] = m.data[5]; r.data[6] = r.data[9];
	r.data[8] = m.data[2]; r.data[9] = m.data[6]; r.data[10]= r.data[10];

	// Inverse Translation
	float tx = m.data[12];
	float ty = m.data[13];
	float tz = m.data[14];

	r.data[12] = -(r.data[0]*tx + r.data[4]*ty + r.data[8]*tz);
	r.data[13] = -(r.data[1]*tx + r.data[5]*ty + r.data[9]*tz);
	r.data[14] = -(r.data[2]*tx + r.data[6]*ty + r.data[10]*tz);

	// letzte Zeile
	r.data[15] = 1.0f;

	return r;
    }

};

template <typename T>
T min (T a, T b) {
    return (a < b) ? a : b;
}



void add_inplace(Vec3& a, const Vec3& b);
Vec3 operator+(const Vec3& a, const Vec3& b);
Vec3 operator-(const Vec3& a, const Vec3& b);
Vec3 operator*(const Vec3& v, float t);
Vec3 operator/(const Vec3& v, float t);
std::ostream& operator<<(std::ostream& outstr, const Vec3 v);
std::string operator<<(const Vec3 v, const char* str);
Vec3 normalize(const Vec3& v);
float dot(Vec3 a, Vec3 b);
Vec3 lerp(Vec3 start, Vec3 end, float t);
float lerpf(float start, float end, float t);
bool float_eq(float a, float b, float eps = 0.0000001f);
float clamp(float t, float low, float high);

Vec3 cross(const Vec3& a, const Vec3& b) {
    return { a.y * b.z - a.z * b.y, 
	     a.z * b.x - a.x * b.z,
	     a.x * b.y - a.y * b.x };
}

template <typename T>
bool in_range(T start, T end, T value) {
    return value < end && value >= start;
}



#ifdef GMATH_IMPLEMENTATION

void add_inplace(Vec3& a, const Vec3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}

Vec3 operator+(const Vec3& a, const Vec3& b) {
    Vec3 v;
    v.x = a.x + b.x;
    v.y = a.y + b.y;
    v.z = a.z + b.z;
    return v;
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    Vec3 v;
    v.x = a.x - b.x;
    v.y = a.y - b.y;
    v.z = a.z - b.z;
    return v;
}

Vec3 operator*(const Vec3& v, float t) {
    Vec3 res;
    res.x = v.x * t;
    res.y = v.y * t;
    res.z = v.z * t;
    return res;
}

Vec3 operator/(const Vec3& v, float t) {
    Vec3 res;
    res.x = v.x / t;
    res.y = v.y / t;
    res.z = v.z / t;
    return res;
}

Mat4 operator*(const Mat4& a, const Mat4& b) {

    Mat4 res;

    for (int row = 0; row < a.size; ++row) {
	int start = row * a.size;

	for (int col = 0; col < a.size; ++col) {
	    res.data[col + row * a.size] = a.data[start] * b.data[col];

	    for (int i = 1; i < a.size; ++i) {
		res.data[col + row * a.size] += a.data[start + i] * b.data[col + a.size * i];
	    }
	}
	
    }
    return res;
}

std::ostream& operator<<(std::ostream& outstr, const Vec3 v) {
	outstr << v.x << " " << v.y << " " << v.z; 
	return outstr;
}

std::string operator<<(const Vec3 v, const char* str) {
	return str;
}

Vec3 normalize(const Vec3& v) {
	return v / v.length();
}

float dot(Vec3 a, Vec3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

float clamp(float t, float low, float high) {
    if (t < low) t = low;
    else if (t > high) t = high;
    return t;
}

Vec3 lerp(Vec3 start, Vec3 end, float t) {
	assert(t >= 0.f && t <= 1.f);
	return start * (1.f - t) + end * t; 
}

float lerpf(float start, float end, float t) {
	assert(t >= 0.f && t <= 1.f);
	//t = clamp(t, 0.f, 1.f);
	return start * (1.f - t) + end * t; 
}

void Vec4::multiply(const Mat4& m) {
    float data[4] = {0};
    int size = 4;
    
    for (int i = 0; i < size; ++i) {
	int start = i * size;
	data[i] = m.data[start] * x + m.data[start + 1] * y + m.data[start + 2] * z + m.data[start + 3] * w;
    }
    x = data[0];
    y = data[1];
    z = data[2];
    w = data[3];
}
void Vec3::multiply(const Mat4& m) {
    float data[3] = {0};
    int size = 3;
    float w = 1.f;

    //std::println("before = {}", to_str());

    for (int i = 0; i < size; ++i) {
	int start_mat = i * m.size;
	data[i] = m.data[start_mat] * x + m.data[start_mat + 1] * y + m.data[start_mat + 2] * z + m.data[start_mat + 3] * w; 
    }
    x = data[0];
    y = data[1];
    z = data[2];

    //std::println("after = {}", to_str());
}

void Vec3::multiply(const Mat3& m) {
    float data[3] = {0};
    int size = 3;

    for (int i = 0; i < size; ++i) {
	int start = i * size;
	data[i] = m.data[start] * x + m.data[start + 1] * y + m.data[start + 2] * z;
    }
    x = data[0];
    y = data[1];
    z = data[2];
}

bool float_eq(float a, float b, float eps) {
    return std::abs(b - a) < eps;
}

#endif //GMATH_IMPLEMENTATION

} // gmath

#endif //GMATH_HPP
