Particle[] planets = new Particle[100];
float[] randomX = new float[100];
float[] randomY = new float[100];
float gravConstScale = 0.6;

class Particle {
  PVector displacement;
  PVector velocity;
  PVector acceleration;
  
  float mass;
  
  Particle (float mass, float x, float y) {
      this.mass = mass;
      this.displacement = new PVector(x, y);
      this.velocity = new PVector(0, 0);
      this.acceleration = new PVector(0, 0);
  }
  
  void applyForce(PVector force) {
    PVector f = force.get();
    f.div(mass);
    acceleration.add(f);
  }
  
  PVector attract(Particle object) {
    
    PVector force = PVector.sub(this.displacement, object.displacement);
    float distance = force.mag();
    distance = constrain(distance, 5.0, 45.0);
    force.normalize();
    float mag = ( gravConstScale * this.mass * object.mass) / (distance * distance);
    
    
    force.mult(mag);
    return force;
    
  }
  
  void checkEdges() {
      if (displacement.x > width) {
        displacement.x = width;
        velocity.x *= -1;
      } else if (displacement.x < 0) {
        velocity.x *= -1;
        displacement.x = 0;
      }
 
      if (displacement.y > height) {
        velocity.y *= -1;
        displacement.y = height;
      } else if (displacement.y < 0) {
        velocity.y *= -1;
        displacement.y = 0;  
      }
  }
  
  void update() {
      velocity.add(acceleration);
      displacement.add(velocity);
      acceleration.mult(0);  
  }
  
  void display() {
    stroke(0);
    fill(175);
    ellipse(displacement.x, displacement.y,mass*10,mass*10);
  }
}

 
void setup() {
  size(400, 400);
 
  
  for (int i = 0; i < planets.length; i++) {
      planets[i] = new Particle(random(0.1, 2), random(width), random(height));  
  }
  
  for (int i = 0; i < randomX.length; i++) {
      randomX[i] = random(width);
  }
  
  for (int i = 0; i < randomY.length; i++) {
      randomY[i] = random(height);
  }
  
  
 
}
 
void draw() {
  background(0);
  
  for (int i = 0; i < 100; i++) {
      fill(255);
      rect(randomX[i], randomY[i], 3, 3);  
  }
     
  for (int i = 0; i < planets.length; i++) {
      for(int j = 0; j < planets.length; j++) {
          if (i != j) {
              PVector gravForce = planets[j].attract(planets[i]);
              planets[i].applyForce(gravForce);
          }
      }
      //planets[i].checkEdges();
      planets[i].update();
      planets[i].display();
  } 
}
