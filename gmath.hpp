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

    std::string to_str() {
	return std::string("Vec3: ") + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z);
    }

    void multiply(const Mat4& m);
    void multiply(const Mat3& m);

    void project() {
	if (z == 0.f) return;
	
	x /= z;
	y /= z;
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
	t = clamp(t, 0.f, 1.f);
	//assert(t >= 0.f && t <= 1.f);
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

    for (int i = 0; i < size; ++i) {
	int start_mat = i * m.size;
	data[i] = m.data[start_mat] * x + m.data[start_mat + 1] * y + m.data[start_mat + 2] * z + m.data[start_mat + 3] * w; 
    }
    x = data[0];
    y = data[1];
    z = data[2];
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
