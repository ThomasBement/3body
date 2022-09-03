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

std::ostream &operator<<(std::ostream &os, vec3 const &m) { 
    return os << "(" << m.x << "," << m.y << "," << m.z << ")";
}

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
    
    void set_matrix(vec3 offset) {
      glm::mat4 model = glm::translate(glm::mat4(1.0), glm::vec3(offset.x, offset.y, offset.z));

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
      set_matrix(vec3(0.0, 0.0, 0.0)); // this will be called many times in the future
    }
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
    size_t getElementsSize() {
        return elements.size() * sizeof(GLuint);
    }
    GLuint *getElements() {
        return elements.data();
    }
    size_t getVerticesSize() {
        return vertices.size() * sizeof(Vert);
    }
    Vert *getVertices() {
        return vertices.data();
    }
    
    void add_cubes(int n) {
        double unl = 1.0/n; // unit length
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
    
    vec3 get_sphere_point(double phi, double theta) {
      return vec3(cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta));
    }
    
    void add_sphere(vec3 center, double r, int N_phi, int N_theta) {
      for (int i = 0; i < N_phi; i++) {
        for (int j = 0; j < N_theta; j++) {
          double phi0 = 2 * M_PI * i / (double)N_phi;
          double phi1 = 2 * M_PI * (i+1) / (double)N_phi;
          double theta0 = M_PI * (j - N_theta/2) / (double)N_theta;
          double theta1 = M_PI * (j + 1 - N_theta/2) / (double)N_theta;
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

vec3 mass_center(vector<planet> &planets) {
  double total_mass = 0.0;
  vec3 sum_p_mass = vec3(0.0, 0.0, 0.0);
  for (int i = 0; i < planets.size(); i++) {
    sum_p_mass = sum_p_mass.plus(planets[i].p.scale(planets[i].mass));
    total_mass += planets[i].mass;
  }
  return sum_p_mass.scale(1/total_mass);
}

double planet_rad(double mass, double dens = 5520) {
  return cbrt((3 * mass) / (4 * M_PI * dens));
}

double rand_dub(double max) {
  return 2 * max * ((double)(rand()) / (double)(RAND_MAX)) - max;
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

void random(vector<planet>& planets, int bodies) {
  for (int i = 0; i < bodies; i++) {
    planet body(0.75 + rand_dub(0.25), vec3(rand_dub(1), rand_dub(1), rand_dub(1)), vec3(rand_dub(0.5), rand_dub(0.5), rand_dub(0.5)));
    planets.push_back(body);
  }
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
    
    srand((unsigned int)time(nullptr));

    // -------------------------------- //
    // PHYSICS CONSTANTS
    // -------------------------------- //
    const double G = 1;                 // Gravitational constant: [m.m.m/kg.s.s] Ours -> 6.6743e-11 
    const double scale = 6e-1;          // Scale factor for visualization 
    
    // -------------------------------- //
    // RENDER CONSTANTS
    // -------------------------------- //
    vec3 orgin_prev = vec3(0.0, 0.0, 0.0);   // For maintaining CG
    vec3 orgin = vec3(0.0, 0.0, 0.0);        // For maintaining CG

    // -------------------------------- //
    // TIME CONSTANTS
    // -------------------------------- //
    double t = 0.0;                     // Physics: Time in while loop [Units of time step]
    const double dt = 1e-4;             // Physics: Time step [s?]
    const int phys_loop = 50;          // Physics: Physics updates per frame [#]
    const double FPS = 60;             // Render: Frames per seccond [FPS]
    const double frame_delay = 1/FPS;   // Render: Frame time period [s]
    glfwSetTime(0.0);                   // Render: Set GLFW timer to zero
    double t0 = 0.0;                    // Render: Previous frame draw time [s]
    double t1 = glfwGetTime();          // Render: Current time [s]

    // Define planet initial conditions and push to vector
    vector<planet> planets = {};
    random(planets, 5);

    // -------------------------------- //
    // WINDOW LOOP
    // -------------------------------- //
    while (glfwWindowShouldClose(mWindow) == false) {
        // -------------------------------- //
        // ESCAPE CONDITIONS
        // -------------------------------- //
        if (glfwGetKey(mWindow, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
          glfwSetWindowShouldClose(mWindow, true);
        }
        
        // -------------------------------- //
        // FRAME UPDATE
        // -------------------------------- //
        t1 = glfwGetTime();
        if ((t1 - t0) >= frame_delay) {
          // -------------------------------- //
          // PHYSICS UPDATE
          // -------------------------------- //
          for (int i = 0; i < phys_loop; i++) {
              phys_up(dt, G, planets);
          }

          // Create world object and add all geometry
          World world;
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

          // Adjust window to follow C.G.
          d.set_matrix(mass_center(planets).scale(-scale));
          
          // Shaders and textures to buffer
          glDrawElements(GL_TRIANGLES, world.getElementsCount(), GL_UNSIGNED_INT, 0);

          // Flip Buffers and Draw
          glfwSwapBuffers(mWindow);
          glfwPollEvents();

          // Update time
          t0 = t1;
        }
    }
    glfwTerminate();
    return EXIT_SUCCESS;
}
