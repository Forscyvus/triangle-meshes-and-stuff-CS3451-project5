// Sample code for starting the meshes project

import processing.opengl.*;

float time = 0;  // keep track of passing of time (for automatic rotation)
boolean rotate_flag = false;       // automatic rotation of model?
boolean randcolors = false;
int seed;
boolean vertexnormals = false;

ArrayList<Vertex> vertices;
ArrayList<Integer> corners;
ArrayList<Integer> opposites;

class Vertex {
  float x, y, z;
  ArrayList<Integer> faces; //holds corner indices (mult by 3)
  ArrayList<Integer> adj; //holds vertex indices
  Vertex(float x, float y, float z){
    this.x = x; this.y = y; this.z = z;
    faces = new ArrayList<Integer>();
    adj = new ArrayList<Integer>();
  }
  void addAdj(int v){
    for (int i = 0; i < adj.size(); i++){
      if (adj.get(i).equals(v)){
        return;
      }
    }
    adj.add(v);
  }
  void addFace(int f){
    faces.add(f);
  }
  Vertex crossProduct(Vertex v){
      return new Vertex(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
  }
  Vertex minus(Vertex v){
      return new Vertex(x-v.x, y-v.y, z-v.z);
  }
  Vertex plus(Vertex v){
      return new Vertex(x+v.x, y+v.y, z+v.z);
  }
  Vertex divide(float d){
    return new Vertex(x/d, y/d, z/d);
  }
  Vertex unit(){
      float length = sqrt(dot(this));
      return new Vertex(x / length, y / length, z / length);
  }
  float dot(Vertex v){
      return (x*v.x + y*v.y + z*v.z);
  }
  boolean equalTo(Vertex v) {
    float tol = .001;
    return ((v.x > (x - tol)) && (v.x < (x + tol)) && (v.y > (y - tol)) && (v.y < (y + tol)) && (v.z > (z - tol)) && (v.z < (z + tol)));
  }
}


// initialize stuff
void setup() {
  size(400, 400, OPENGL);  // must use OPENGL here !!!
  noStroke();     // do not draw the edges of polygons
  vertices = new ArrayList<Vertex>();
  corners = new ArrayList<Integer>();
  opposites = new ArrayList<Integer>();
}

// Draw the scene
void draw() {
  
  resetMatrix();  // set the transformation matrix to the identity (important!)

  background(0);  // clear the screen to black
  
  // set up for perspective projection
  perspective (PI * 0.333, 1.0, 0.01, 1000.0);
  
  // place the camera in the scene (just like gluLookAt())
  camera (0.0, 0.0, 5.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
  
  scale (1.0, -1.0, 1.0);  // change to right-handed coordinate system
  
  // create an ambient light source
  ambientLight(102, 102, 102);
  
  // create two directional light sources
  lightSpecular(204, 204, 204);
  directionalLight(102, 102, 102, -0.7, -0.7, -1);
  directionalLight(152, 152, 152, 0, 0, -1);
  
  pushMatrix();

  fill(200, 200, 200);            // set polygon color to blue
  ambient (200, 200, 200);
  specular(0, 0, 0);
  shininess(1.0);
  
  rotate (time, 1.0, 0.0, 0.0);
  
  // THIS IS WHERE YOU SHOULD DRAW THE MESH
  
  randomSeed(seed);
  for (int j = 0; j < corners.size()/3; j++) {
    Vertex v1 = vertices.get(corners.get(3*j));
    Vertex v2 = vertices.get(corners.get(3*j + 1));
    Vertex v3 = vertices.get(corners.get(3*j + 2));
    beginShape();
    
    if(randcolors){
      fill(random(0,255),random(0,255),random(0,255));
    }
    if (vertexnormals) {
      Vertex norm = getNormal(v1);
      normal(norm.x, norm.y, norm.z);
    }
    vertex(v1.x, v1.y, v1.z);
    if (vertexnormals) {
      Vertex norm = getNormal(v2);
      normal(norm.x, norm.y, norm.z);
    }
    vertex(v2.x, v2.y, v2.z);
    if (vertexnormals) {
      Vertex norm = getNormal(v3);
      normal(norm.x, norm.y, norm.z);
    }
    vertex(v3.x, v3.y, v3.z);
    endShape(CLOSE);
  }
  
  popMatrix();
 
  // maybe step forward in time (for object rotation)
  if (rotate_flag )
    time += 0.02;
}

// handle keyboard input
void keyPressed() {
  if (key == '1') {
    read_mesh ("tetra.ply");
  }
  else if (key == '2') {
    read_mesh ("octa.ply");
  }
  else if (key == '3') {
    read_mesh ("icos.ply");
  }
  else if (key == '4') {
    read_mesh ("star.ply");
  }
  else if (key == '5') {
    read_mesh ("torus.ply");
  }
  else if (key == '6') {
    create_sphere(5,16);                     // create a sphere
  }
  else if (key == ' ') {
    rotate_flag = !rotate_flag;          // rotate the model?
  }
  else if (key == 'q' || key == 'Q') {
    exit();                               // quit the program
  }
  else if (key == 'r' || key == 'R') {
    randcolors = true;
    seed = (int) random(0,1000);
  }
  else if (key == 'w' || key == 'W') {
    randcolors = false;
  }
  else if (key == 'n' || key == 'N') {
    vertexnormals = !vertexnormals;
  }
  else if (key == '0') {
    for (int i = 0; i < vertices.size(); i++){
      Vertex n = getNormal(vertices.get(i));
      println("Normal " + i + ": " + n.x + ' ' + n.y + ' ' + n.z);
    }
  }
  else if (key == 'd' || key == 'D') {
    dual();
  }
}

Vertex getNormal(Vertex v){
  
  Vertex result = new Vertex(0.0,0.0,0.0);
  
  for (Integer i : v.faces){
    Vertex v1 = vertices.get(corners.get(i*3));
    Vertex v2 = vertices.get(corners.get(i*3 + 1));
    Vertex v3 = vertices.get(corners.get(i*3 + 2));
    Vertex norm = v2.minus(v1).unit().crossProduct(v3.minus(v1).unit());
    
    result = result.plus(norm);
  }
  
  return result.divide(v.faces.size()).unit();
  //return new Vertex(0,0,1);
}

// Read polygon mesh from .ply file
//
// You should modify this routine to store all of the mesh data
// into a mesh data structure instead of printing it to the screen.
void read_mesh (String filename)
{
  vertices = new ArrayList<Vertex>();
  corners = new ArrayList<Integer>();
  
  int i;
  String[] words;
  
  String lines[] = loadStrings(filename);
  
  words = split (lines[0], " ");
  int num_vertices = int(words[1]);
  println ("number of vertices = " + num_vertices);
  
  words = split (lines[1], " ");
  int num_faces = int(words[1]);
  println ("number of faces = " + num_faces);
  
  // read in the vertices
  for (i = 0; i < num_vertices; i++) {
    words = split (lines[i+2], " ");
    float x = float(words[0]);
    float y = float(words[1]);
    float z = float(words[2]);
    println ("vertex = " + x + " " + y + " " + z);
    vertices.add(new Vertex(x,y,z));
  }
  
  // read in the faces
  for (i = 0; i < num_faces; i++) {
    
    int j = i + num_vertices + 2;
    words = split (lines[j], " ");
    
    int nverts = int(words[0]);
    if (nverts != 3) {
      println ("error: this face is not a triangle.");
      exit();
    }
    
    int index1 = int(words[1]);
    int index2 = int(words[2]);
    int index3 = int(words[3]);
    println ("face = " + index1 + " " + index2 + " " + index3);
    corners.add(index1);
    corners.add(index2);
    corners.add(index3);
//    vertices.get(index1).addFace(i);
//    vertices.get(index2).addFace(i);
//    vertices.get(index3).addFace(i);
  }
  generateAdjacencies();
}

void generateAdjacencies(){
  //update amount of faces per vertex
  for (int i = 0; i < corners.size(); i += 3){
    vertices.get(corners.get(i)).addFace(i/3);
    vertices.get(corners.get(i+1)).addFace(i/3);
    vertices.get(corners.get(i+2)).addFace(i/3);
    
    vertices.get(corners.get(i)).addAdj(i+1);
    vertices.get(corners.get(i)).addAdj(i+2);
    vertices.get(corners.get(i+1)).addAdj(i);
    vertices.get(corners.get(i+1)).addAdj(i+2);
    vertices.get(corners.get(i+2)).addAdj(i);
    vertices.get(corners.get(i+2)).addAdj(i+1);
  }
  
  //do actual corner stuff
  opposites = new ArrayList<Integer>();
  for (int i = 0; i < corners.size(); i++){
    opposites.add(null);
  }
  for (int i = 0; i < corners.size(); i++){
    int in = 0; 
    int ip = 0;
    if (i % 3 == 0){
      in = i + 1; ip = i + 2;
    } else if (i % 3 == 1) {
      in = i + 1; ip = i - 1;
    } else {
      in = i - 2; ip = i - 1;
    }
    for (int j = 0; j < corners.size(); j++){
      int jn = 0; 
      int jp = 0;
      if (j % 3 == 0){
        jn = j + 1; jp = j + 2;
      } else if (j % 3 == 1) {
        jn = j + 1; jp = j - 1;
      } else {
        jn = j - 2; jp = j - 1;
      }
      
      if( (corners.get(in).equals(corners.get(jp))) && (corners.get(jn).equals( corners.get(ip)) )){
        opposites.set(i,j);
        opposites.set(j,i);
      }
    }
  }
  
  println("opposites: " + opposites);
  //println("corners: " + corners);
  
//  ArrayList<Integer> counts = new ArrayList<Integer>();
//  for (int i = 0; i < vertices.size(); i++){
//    counts.add(0);
//  }
//  for (int i = 0; i < corners.size(); i++){
//    int v = corners.get(i);
//    counts.set(v, counts.get(v) + 1);
//  }
//  println("counts: " + counts);
  
}

//lats = latitudes per hemisphere
//longs = longitudes all the way around
void create_sphere(int lats, int longs){

  vertices = new ArrayList<Vertex>();
  corners = new ArrayList<Integer>();

  //vertices
  for (int lat = 1 - lats; lat < lats; lat ++){
    for (int lon = 0; lon < longs; lon ++){
      
      float t = HALF_PI*(lat/ ((float)lats) );
      float s = TWO_PI*(lon/((float)longs));
      float tplus = HALF_PI*((lat+1)/ ((float)lats) );
      float tminus = HALF_PI*((lat-1)/ ((float)lats) );
      float splus = TWO_PI*((lon+1)/((float)longs));
      
      //upper face
      corners.add( getVertexIndex( cos(t) * cos(s), cos(t) * sin(s), sin(t), vertices ) );
      corners.add( getVertexIndex( cos(t) * cos(splus), cos(t) * sin(splus), sin(t), vertices ) );
      corners.add( getVertexIndex( cos(tplus) * cos(s), cos(tplus) * sin(s), sin(tplus), vertices ) );
      
      //lower face
      corners.add( getVertexIndex( cos(t) * cos(s), cos(t) * sin(s), sin(t), vertices ) );
      corners.add( getVertexIndex( cos(tminus) * cos(splus), cos(tminus) * sin(splus), sin(tminus), vertices ) );
      corners.add( getVertexIndex( cos(t) * cos(splus), cos(t) * sin(splus), sin(t), vertices ) );
      
      
      
      }
    }
    
    generateAdjacencies();
    
//    for (int i = 0; i < vertices.size(); i++) {
//      println("Vertex: ");
//      println(vertices.get(i).x + " " + vertices.get(i).y + " " + vertices.get(i).z);
//    }
    
    println("Vertices Size: " );
    println(vertices.size());
    
//    println("Corners: ");
//    println(corners);
//    
//    println("Corners Size: " );
//    println(corners.size()/3); 
  }
  
int getVertexIndex(float x, float y, float z, ArrayList<Vertex> list) {
  for (int i = 0; i < list.size(); i++) {
    Vertex v = list.get(i);
    
//    if (i == 1) {
//      println("xyz, v.xyz: " + x + " " + v.x + " / " + y + " " + v.y + " / " + z + " " + v.z + " / ");
//    }
    
    //if ((v.x > (x - tol)) && (v.x < (x + tol)) && (v.y > (y - tol)) && (v.y < (y + tol)) && (v.z > (z - tol)) && (v.z < (z + tol))){
    //if (v.x == x && v.y == y && v.z == z){
    if (v.equalTo(new Vertex(x,y,z))){
      return i;
    }
  }
  list.add(new Vertex(x,y,z));
  return list.size()-1;  
}
  
void dual(){
  ArrayList<Vertex> newvertices = new ArrayList<Vertex>();
  ArrayList<Integer> newcorners = new ArrayList<Integer>();


  for (int i = 0; i < vertices.size(); i++) {
    Vertex v = vertices.get(i);
    ArrayList<Vertex> cents = new ArrayList<Vertex>();
    int startface = v.faces.get(0);
    int face = startface;
    int failsafe = 0;
    do {
      cents.add(centroid(face));
      int corneroffset = 0;
      for( int j = 0; j < 3; j++){
        if ( corners.get(face*3 + j) == i ){
          corneroffset = j;
        }
      }
      face = opposites.get( (face*3 + ((corneroffset+1)%3)) ) / 3;
      failsafe++;
    }while( face != startface && failsafe < 100);
    
    Vertex newcentroid = centroid(cents);
    for (int j = 0; j < cents.size(); j++) {
      int k = (j + 1)%cents.size();
      newcorners.add( getVertexIndex(newcentroid.x, newcentroid.y, newcentroid.z, newvertices));
      newcorners.add( getVertexIndex(cents.get(j).x, cents.get(j).y, cents.get(j).z, newvertices));
      newcorners.add( getVertexIndex(cents.get(k).x, cents.get(k).y, cents.get(k).z, newvertices));
    }
    
  }
  
  
  
  
  vertices = newvertices;
  corners = newcorners;
  generateAdjacencies();
}

Vertex centroid( ArrayList<Vertex> poly ) {
  
  float x = 0; float y = 0; float z = 0;
  for (int i = 0; i < poly.size(); i++ ) {
    x += poly.get(i).x;
    y += poly.get(i).y;
    z += poly.get(i).z; 
  }
  
  x = x / (float)poly.size();
  y = y / (float)poly.size();
  z = z / (float)poly.size();
  
  return new Vertex(x,y,z);
  
}
Vertex centroid( int t ) {
  ArrayList<Vertex> triangle = new ArrayList<Vertex>();
  triangle.add(vertices.get(corners.get(t*3)));
  triangle.add(vertices.get(corners.get(t*3 + 1)));
  triangle.add(vertices.get(corners.get(t*3 + 2)));
  return centroid(triangle);
}

