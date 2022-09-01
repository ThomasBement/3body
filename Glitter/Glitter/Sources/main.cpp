// Defines
#define _USE_MATH_DEFINES

// Local Headers
#include "glitter.hpp"

// System Headers
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// Math
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cmath>

// Standard Headers
#include <cstdio>
#include <cstdlib>
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout; using std::endl;
using std::string;
#include <vector>
using std::vector;
#include <thread>
using namespace std::this_thread;
#include <chrono>
using namespace std::chrono;

// Images
#include "stb_image.h"

string read_file(const string filename) {
  ifstream f(filename);
  if (f) {
    std::string s;
    f.seekg(0, std::ios::end);
    const int size = f.tellg();
    s.resize(size);
    f.seekg(0, std::ios::beg);
    f.read((char*)s.data(), size);
    f.close();
    return s;
  } else throw(errno);
}

void log_shader_info(GLuint shader) {
  char buffer[512];
  glGetShaderInfoLog(shader, 512, NULL, buffer);
  fprintf(stderr, "%s\n", buffer);
}

void log_program_info(GLuint program) {
  char buffer[512];
  glGetProgramInfoLog(program, 512, NULL, buffer);
  fprintf(stderr, "%s\n", buffer);
}

GLuint load_shader_source(const string filename, GLenum type) {
  string vs_source = read_file(filename);
  const GLchar *data = vs_source.data();
  const GLint length = vs_source.size();
  GLuint shader = glCreateShader(type);
  glShaderSource(shader, 1, &data, &length);
  glCompileShader(shader);
  fprintf(stderr, "Compilation log for shader %s:\n", filename.c_str());
  log_shader_info(shader);
  return shader;
}

class D {
  private:
    GLuint vbo;
    GLuint vertexShader;
    GLuint fragmentShader;
    GLuint shaderProgram;
    GLint posAttrib;
    GLint colAttrib;
    GLuint vao;
    GLint uniColor;
    GLuint ebo;
    GLuint tex;
    GLint texAttrib;
    
  public:
    void vertexData(GLsizeiptr size, const GLvoid *vertices) {
      glBufferData(GL_ARRAY_BUFFER, size, vertices, GL_DYNAMIC_DRAW);
    }
    
    void elementData(GLsizeiptr size, const GLvoid *elements) {
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, elements, GL_DYNAMIC_DRAW);
    }
    
    void load_image_texture(string filename) {
      int width; int height; int channels;
      uint8_t *pixels = stbi_load(filename.c_str(), &width, &height, &channels, 0);
      if (stbi_failure_reason()) {
        fprintf(stderr, "STB Image Lib Error: %s\n", stbi_failure_reason());
      } else if (channels != 3) {
        fprintf(stderr, "Error: %d channels in image %s, expected 3.\n",
          channels, filename.c_str());
      }
      fprintf(stderr, "hello world!\n");
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
      stbi_image_free(pixels);
    }
    
    void projection_matrix_setup() {
      float aspect_ratio = (float)mWidth / (float)mHeight;
      glm::mat4 proj = glm::perspective(glm::radians(45.0f), aspect_ratio, 0.1f, 10.0f);
      GLint uniProj = glGetUniformLocation(shaderProgram, "proj");
      glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));
      // want to allow depth  buffer and depth comparison:
      glEnable(GL_DEPTH_TEST);
    }
    
    void set_matrix(float t) {
      glm::mat4 model = glm::mat4(1.0f);
      model = glm::rotate(model, glm::radians(t), glm::vec3(0.0f, 0.0f, 1.0f));
      GLint uniModel = glGetUniformLocation(shaderProgram, "model");
      glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
      glm::mat4 view = glm::lookAt(
          glm::vec3(1.2f, 1.2f, 1.2f),
          glm::vec3(0.0f, 0.0f, 0.0f),
          glm::vec3(0.0f, 0.0f, 1.0f)
      );
      GLint uniView = glGetUniformLocation(shaderProgram, "view");
      glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
    }
    
    void init_vertex_attributes() {
      posAttrib = glGetAttribLocation(shaderProgram, "position");
      glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), 0);
      glEnableVertexAttribArray(posAttrib);
      colAttrib = glGetAttribLocation(shaderProgram, "color");
      glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*)(3*sizeof(float)));
      glEnableVertexAttribArray(colAttrib);
      texAttrib = glGetAttribLocation(shaderProgram, "texcoord");
      glVertexAttribPointer(texAttrib, 2, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*)(6*sizeof(float)));
      glEnableVertexAttribArray(texAttrib);
    }
    
    D() {
      // vertex array object
      glGenVertexArrays(1, &vao);
      glBindVertexArray(vao);
      // vertex buffer
      glGenBuffers(1, &vbo);
      glBindBuffer(GL_ARRAY_BUFFER, vbo);
      // element buffer
      glGenBuffers(1, &ebo);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
      // textures
      glGenTextures(1, &tex);
      glBindTexture(GL_TEXTURE_2D, tex);
      load_image_texture("assets/textures.jpeg");
      fprintf(stderr, "Loaded texture pixels. Error code is: %d \n", glGetError());
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      fprintf(stderr, "Set out of bounds to repeat. Error code is: %d \n", glGetError());
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      fprintf(stderr, "Loaded texture params. Error code is: %d \n", glGetError());
      // shaders
      vertexShader   = load_shader_source("vertex.glsl", GL_VERTEX_SHADER);
      fragmentShader = load_shader_source("fragment.glsl", GL_FRAGMENT_SHADER);
      shaderProgram = glCreateProgram();
      glAttachShader(shaderProgram, vertexShader);
      glAttachShader(shaderProgram, fragmentShader);
      glBindFragDataLocation(shaderProgram, 0, "outColor");
      glLinkProgram(shaderProgram);
      glUseProgram(shaderProgram);
      // vertex attributes
      init_vertex_attributes();
      // uniforms
      uniColor = glGetUniformLocation(shaderProgram, "triangleColor");
      glUniform3f(uniColor, 0.5f, 0.0f, 1.0f);
      projection_matrix_setup(); // we only call this once
      set_matrix(0.0); // this will be called many times in the future
    }
};

class vec3 {
  public:
    double x;
    double y;
    double z;
    
    vec3(double xx, double yy, double zz)
      : x(xx), y(yy), z(zz) {}
    
    vec3 plus(vec3 v) {
        return vec3(x+v.x, y+v.y, z+v.z);
    }

    vec3 sub(vec3 v) {
        return vec3(x-v.x, y-v.y, z-v.z);
    }

    vec3 scale(double a) {
        return vec3(a*x, a*y, a*z);
    }

    double abs() {
        return sqrt(x*x + y*y + z*z);
    }
};

class planet {       
    public:            
        double mass;
        double scale;
        vec3 p;
        vec3 v;
        vec3 a;
        vec3 a0;
        planet(double mm, vec3 pp, vec3 vv, double ss = 1) 
            : mass(mm), scale(ss), p(pp), v(vv), a(0.0, 0.0, 0.0), a0(0.0, 0.0, 0.0) {}
    
    void p_up(double dt) {
        this->p = this->p.plus(this->v.scale(dt));
    }

    void v_up(double dt) {
        this->v = this->v.plus(this->a.scale(dt));
    }
};

class Vert {
  public:
    GLfloat x;
    GLfloat y;
    GLfloat z;
    GLfloat r = 1.0;
    GLfloat g = 1.0;
    GLfloat b = 1.0;
    GLfloat u;
    GLfloat v;
    
    Vert(float xx, float yy, float zz, float uu, float vv)
      : x(xx), y(yy), z(zz), u(uu), v(vv) {}
      
    Vert(vec3 pos, float uu, float vv)
      : x(pos.x), y(pos.y), z(pos.z), u(uu), v(vv) {}
};

class World {
  private:
    vector<Vert> vertices;
    vector<GLuint> elements;
    
    void add_face(vec3 base, vec3 du, vec3 dv) {
        int vertexLoc = vertices.size();
        vertices.push_back(Vert(base, 3., 0.));
        vertices.push_back(Vert(base.plus(du), 4., 0.));
        vertices.push_back(Vert(base.plus(dv), 3., 1.));
        vertices.push_back(Vert(base.plus(du).plus(dv), 4., 1.));
        elements.push_back(vertexLoc);
        elements.push_back(vertexLoc+1);
        elements.push_back(vertexLoc+2);
        elements.push_back(vertexLoc+1);
        elements.push_back(vertexLoc+2);
        elements.push_back(vertexLoc+3);
    }
    
    void add_quad(vec3 a, vec3 b, vec3 c, vec3 d) {
      int vertexLoc = vertices.size();
      vertices.push_back(Vert(a, 3., 0.));
      vertices.push_back(Vert(b, 3., 1.));
      vertices.push_back(Vert(c, 4., 0.));
      vertices.push_back(Vert(d, 4., 1.));
      elements.push_back(vertexLoc);
      elements.push_back(vertexLoc+1);
      elements.push_back(vertexLoc+2);
      elements.push_back(vertexLoc+1);
      elements.push_back(vertexLoc+2);
      elements.push_back(vertexLoc+3);
    }
    
  public:    
    int getElementsCount() {
        return elements.size();
    }
    int getElementsSize() {
        return elements.size() * sizeof(GLuint);
    }
    GLuint *getElements() {
        return elements.data();
    }
    int getVerticesSize() {
        return vertices.size() * sizeof(Vert);
    }
    Vert *getVertices() {
        return vertices.data();
    }
    
    void add_cubes(int n) {
        float unl = 1.0/n; // unit length
        vec3 dx(unl, 0.0, 0.0);
        vec3 dy(0.0, unl, 0.0);
        vec3 dz(0.0, 0.0, unl);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i + j) % 2 == 0) {
                    add_face(vec3(i*unl, j*unl, 0.0), dx, dy);
                    add_face(vec3(i*unl, j*unl, 0.0), dy, dz);
                    add_face(vec3(i*unl, j*unl, 0.0), dz, dx);
                    add_face(vec3(i*unl, j*unl, unl), dx, dy);
                    add_face(vec3((i+1)*unl, j*unl, 0.0), dy, dz);
                    add_face(vec3(i*unl, (i+1)*unl, 0.0), dz, dx);
                }
            }
        }
    }
    
    vec3 get_sphere_point(float phi, float theta) {
      return vec3(cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta));
    }
    
    void add_sphere(vec3 center, float r, int N_phi, int N_theta) {
      for (int i = 0; i < N_phi; i++) {
        for (int j = 0; j < N_theta; j++) {
          float phi0 = 2 * M_PI * i / (float)N_phi;
          float phi1 = 2 * M_PI * (i+1) / (float)N_phi;
          float theta0 = M_PI * (j - N_theta/2) / (float)N_theta;
          float theta1 = M_PI * (j + 1 - N_theta/2) / (float)N_theta;
          vec3 a = center.plus(get_sphere_point(phi0, theta0).scale(r));
          vec3 b = center.plus(get_sphere_point(phi0, theta1).scale(r));
          vec3 c = center.plus(get_sphere_point(phi1, theta0).scale(r));
          vec3 d = center.plus(get_sphere_point(phi1, theta1).scale(r));
          add_quad(a, b, c, d);
        }
      }
    }
  
    World()
      : vertices(), elements() {}
};

void phys_up(double dt, double G, vector<planet> &planets) {
    // Position update: Vec3d new_pos = pos + vel * dt + acc * (dt * dt * 0.5);
    for (int i = 0; i < planets.size(); i++) {
        planets[i].p = planets[i].p.plus(planets[i].v.scale(dt)).plus(planets[i].a.scale(dt * dt * 0.5));
        //planets[i].p = planets[i].p.plus(planets[i].v.plus(planets[i].a.scale(0.5 * dt)).scale(dt));
    }
    // Save previous acceleration
    for (int i = 0; i < planets.size(); i++) {
      planets[i].a0 = planets[i].a;
      planets[i].a = vec3(0.0, 0.0, 0.0);
    }
    // Compute acceleration
    for (int i = 0; i < planets.size(); i++) {
        planets[i].a0 = planets[i].a;
        for (int j = 0; j < i; j++) {
            vec3 displacement = planets[i].p.sub(planets[j].p);
            double distance = displacement.abs();
            vec3 G_inverse_square_ij = displacement.scale(G / (distance * distance * distance));
            planets[i].a = planets[i].a.plus(G_inverse_square_ij.scale(-planets[j].mass));
            planets[j].a = planets[j].a.plus(G_inverse_square_ij.scale(planets[i].mass));
        }
    }
    // Velocity update: Vec3d new_vel = vel + (acc + new_acc) * (dt * 0.5);
    for (int i = 0; i < planets.size(); i++) {
        planets[i].v = planets[i].v.plus((planets[i].a0.plus(planets[i].a)).scale(dt * 0.5));
    }
}

double planet_rad(double mass, double dens = 5520) {
    return cbrt((3 * mass) / (4 * M_PI * dens));
}

void lagrange_system(vector<planet> &planets) {
    const double rad_sep = 1;
    const double v_scale = 5e-1;
    
    vector<vec3> p_init = {};
    vector<vec3> v_init = {};
    
    for (int i = 0; i < 3; i++) {
        double theta = 2 * i * M_PI / 3;
        p_init.push_back(vec3(rad_sep * cos(theta), rad_sep * sin(theta), 0));
        v_init.push_back(vec3(v_scale * rad_sep * sin(theta), -v_scale * rad_sep * cos(theta), 0));
    }

    planet sun_1(1, p_init[0], v_init[0]);
    planet sun_2(1, p_init[1], v_init[1]);
    planet sun_3(1, p_init[2], v_init[2]);
    planet planet_1(1.0/1000, vec3(5e-1, -5e-1, 2e-1), vec3(0, 5e-1, 0), 5e0);

    planets.push_back(sun_1);
    planets.push_back(sun_2);
    planets.push_back(sun_3);
    planets.push_back(planet_1);
}

void figure_eight(vector<planet>& planets) {
    planet sun_1(1, vec3(-0.97000436, 0.24308753, 0), vec3(0.466203685, 0.43236573, 0));
    planet sun_2(1, vec3(0.97000436, -0.24308753, 0), vec3(0.466203685, 0.43236573, 0));
    planet sun_3(1, vec3(0, 0, 0), vec3(-2*0.466203685, -2*0.43236573, 0));

    planets.push_back(sun_1);
    planets.push_back(sun_2);
    planets.push_back(sun_3);
}

int main() {

    // Load GLFW and Create a Window
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    auto mWindow = glfwCreateWindow(mWidth, mHeight, "OpenGL", nullptr, nullptr);

    // Check for Valid Context
    if (mWindow == nullptr) {
        fprintf(stderr, "Failed to Create OpenGL Context");
        return EXIT_FAILURE;
    }

    // Create Context and Load OpenGL Functions
    glfwMakeContextCurrent(mWindow);
    gladLoadGL();
    fprintf(stderr, "OpenGL %s\n", glGetString(GL_VERSION));

    // Create and initialize the drawing object:
    D d;
    fprintf(stderr, "Loaded vertex data. Error code is: %d \n", glGetError());
    const double G = 1; // 6.6743e-11 [m.m.m/kg.s.s]
    const double scale = 6e-1;
    
    // Physics time
    double t = 0.0;
    const double dt = 1e-4;
    const int phys_loop = 1e2;
    
    // Display time
    const int FPS = 120;
    const int nano_delay = 1e9 / FPS;

    // Define planet initial conditions and push to vector
    vector<planet> planets = {};
    lagrange_system(planets);

    // Rendering Loop
    while (glfwWindowShouldClose(mWindow) == false) {
        if (glfwGetKey(mWindow, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(mWindow, true);
        
        // Create the first world objects based on planets vector
        World world; // New frame new world
        for (int i = 0; i < planets.size(); i++) {
          world.add_sphere(planets[i].p.scale(scale), scale*planet_rad(planets[i].mass)*planets[i].scale, 16, 8);
        }

        // Get verticies and elements for all objects in world
        d.vertexData(world.getVerticesSize(), world.getVertices());
        d.elementData(world.getElementsSize(), world.getElements());

        // Background Fill Color
        glClearColor(0, 0, 0, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        t += dt;
        
        for (int i = 0; i < phys_loop; i++) {
            phys_up(dt, G, planets);
        }

        //d.set_matrix(t); View rotation
        
        // Shaders and textures to buffer
        glDrawElements(GL_TRIANGLES, world.getElementsCount(), GL_UNSIGNED_INT, 0);

        // Flip Buffers and Draw
        glfwSwapBuffers(mWindow);
        glfwPollEvents();

        // Apply FPS delay
        sleep_for(nanoseconds(nano_delay));
    }
    glfwTerminate();
    return EXIT_SUCCESS;
}
