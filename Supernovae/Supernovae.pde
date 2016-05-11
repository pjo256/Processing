/* Maps type Ib, Ic, IIa supernovae from the X survey to a model of the celestial sphere
 * Each event is mapped to the surface of the sphere using its right ascension and declination
 * A vector perpendicular to this point, with a length scaled by the apparent brightness of the supernova
 * The vector is drawn using a color gradient to present spectral types
 * 
*/

import java.text.DateFormat;
//import java.awt.event.*;
/* Global variables */

// Latitude / Longitude on celestial sphere (in radians)
float lat, lon = 0;
float scale_val = 1.0;
float delta_lat = 0.0;
float delta_lon = 0.0;
PFont data;

double mean_r_band;
double mean_v_band;
double mean_i_band;

double[] v_stats; //Store minimum [0] and maximum [1] differences in v - r bands and r - i bands
double[] i_stats;
double[] r_stats;

boolean show_data;

PImage bg;
PImage texmap;

int sDetail = 35;  // Sphere detail setting
float rotationX = 0;
float rotationY = 0;
float velocityX = 0;
float velocityY = 0;
float pushBack = 0;

float[] cx, cz, sphereX, sphereY, sphereZ;
float sinLUT[];
float cosLUT[];
float SINCOS_PRECISION = 0.5;
int SINCOS_LENGTH = int(360.0 / SINCOS_PRECISION);



//double mean_r_i_band;

//Number of fields obtained from SNinfo.txt
int numInitProperties = 14;
SNova events[];
SNova binary_events[]; //Supernovae from binary systems

float celestial_radius = 120;
float globeRadius = celestial_radius;

float mouseZoom;

public class SNova {
  

  String name;
  String type;
  String hostGalaxy;
  
  String formattedRA;
  String formattedDec;
  
  String formattedDate;
  
  float RA;
  float dec;
  
  float julianDate;
  float v_band_mag; //V filter
  float r_band_mag; //R filter
  float i_band_mag; //I filter
  
  int[] rgb_colors;
  
  PVector vect;
  float xPos;
  float yPos;
  float zPos;
  
  //Diameter of sphere containing info
  float radius = 1;
  //boolean captured = false;
  
  public SNova(String name, String type, String hostGalaxy, String formattedRA, String formattedDec, float RA, float dec) {
    this.name = name;
    this.type = type;
    this.hostGalaxy = hostGalaxy;
    this.formattedRA = formattedRA;
    this.formattedDec = formattedDec;
    this.RA = RA;
    this.dec = dec;
    rgb_colors = new int[3];
    
  }
  
  
  
  
  void renderText() {
    noStroke();
    textAlign(LEFT);
    String info = formatText();
    float w = textWidth(info);
    fill(255);
    text(info, 0, 0); 
  }
  
  
  void showMouseOver(PVector data_line) {
    float xPos = screenX(data_line.x, data_line.y, data_line.z);
    float yPos = screenY(data_line.x, data_line.y, data_line.z) ;
    float zPos = screenZ(data_line.x, data_line.y, data_line.z);
    
    //Only show text when cursor is inside ellipse
    float diff_x = mouseX - xPos;
    float diff_y = mouseY - yPos;
    
    
    pushMatrix();
    translate(data_line.x, data_line.y, data_line.z);
    
    if((sq(diff_x) + sq(diff_y) < sq(radius) + 10) && !show_data){
        radius = 4;
        show_data = true;
    } else {
        show_data = false;
        radius = 1;  
    }
    ellipse(0, 0, 2 * radius, 2 * radius);
    if (show_data) {
      
      renderText();  
    }
    popMatrix();
  
  }
  
  
  
  private String formatText() { 
    String info = "Name: " + name + "\n" + "Type: " + type + "\n";
    
    if (this.hostGalaxy != " ") {
      info += "Host Galaxy: " + this.hostGalaxy + "\n"; 
    }
    
    
    info += "RA: " + this.formattedRA + "\n" + "Dec: " + this.formattedDec + "\n" + "Discovery Time (Jul.): " + this.julianDate + "Discovery Date (Greg.): " + this.formattedDate;
    return info;
  }
  
  private void setCartesian(float RA, float dec, float radius) {
    
    float polar_dec = dec + HALF_PI;                     // shift dec by 90 degrees to get polar angle
    
     
    this.xPos = -radius * sin(polar_dec) * cos(RA);     // -1 to correct orientation of Processing's P3D coordinates
    this.zPos =  radius * sin(polar_dec) * sin(RA);
    this.yPos =  radius * cos(polar_dec);
   
    
    this.vect = new PVector(this.xPos, this.yPos, this.zPos);
    
  }  
  
  private void setDateMagnitudes(float julianDate, float v_band_mag, float r_band_mag, float i_band_mag) {
    this.julianDate = julianDate;
    this.v_band_mag = v_band_mag;
    this.r_band_mag = r_band_mag;
    this.i_band_mag = i_band_mag;  
    
  }
  
  /*
  private void setGregorianDate(float julianDate) {
    int JGREG= 15 + 31*(10+12*1582);
    int jalpha, ja, jb, jc, jd, je, year, month, day;
    
    double julian = julianDate + (0.5 / 86400.0);
    ja = (int) julian;
    if (ja>= JGREG) {
      jalpha = (int) (((ja - 1867216) - 0.25) / 36524.25);
      ja = ja + 1 + jalpha - jalpha / 4;
    }

    jb = ja + 1524;
    jc = (int) (6680.0 + ((jb - 2439870) - 122.1) / 365.25);
    jd = 365 * jc + jc / 4;
    je = (int) ((jb - jd) / 30.6001);
    day = jb - jd - (int) (30.6001 * je);
    month = je - 1;
    
    if (month > 12) { 
      month = month - 12;
    }
    
    year = jc - 4715;
    if (month > 2) { year--; }
    if (year <= 0) { year--; }

    this.formattedDate = day + "/" + month + "/" + year;
  }
  */
  
  private void computeColor() {
    
    
    //Maps v_band to red, r_band to green, and i_band to blue
    //Since inputs represent apparent magnitudes, v_band > i_band implies that an event is less bright in the UV end of the spectra
    //Compute v - r, r - i. 
    //If no r-band, set r to average of blue and red.
    //rgb_colors = new int[3];
    
    
    float v_band_mag = this.v_band_mag;
    float r_band_mag = this.r_band_mag;
    float i_band_mag = this.i_band_mag;
    
    
    if (v_band_mag == 0 && i_band_mag == 0) {
       
      v_band_mag = (int) mean_v_band;
      i_band_mag = (int) mean_i_band;
      
    }  else if (v_band_mag == 0) {
      
      v_band_mag = (int) mean_v_band;
    } else if (i_band_mag == 0) {
      i_band_mag = (int) mean_i_band;
    } 

      //Redness
      
      //double max = v_r_diff[1];
      //double min = v_r_diff[0];
      
      double max = v_stats[1];
      double min = v_stats[0];
      rgb_colors[0] = 255 - (int) ((v_band_mag - min)*255/(max-min));
      //Blueness
      max = i_stats[1];
      min = i_stats[0];
       
      rgb_colors[2] = 255 - (int) ((i_band_mag - min)*255/(max-min));
      
       //Fixed green
      max = r_stats[1];
      min = r_stats[0];
      rgb_colors[1] =  (int) ((126 - ((mean_r_band - min) * 126  / (max - min))) - (Math.abs(rgb_colors[0] - rgb_colors[2]) * (mean_r_band - min)));
      
    
  }
}

/* Helpers */

public float getRA(String input) {
  //Input: hh:mm:ss
  //Output: rads
  
  //Remove semi-colons
  
  String normalized[] = new String[3];
  normalized = input.split(":");
  
  
  double degrees = Double.parseDouble(normalized[0]) * 15.0;
  
  double minToDeg = Double.parseDouble(normalized[1]) * (15.0 / 60.0);
  double secToDeg = Double.parseDouble(normalized[2]) * (15.0 / 3600.0);
 
  
  degrees += (minToDeg + secToDeg);
  
  double rads = ((Math.PI) / (180.0)) * degrees;
  
  return (float)rads;
    
}

public float getDec(String input) {
  //Input: +/- dd:mm:ss
  //Output: rads
  
  String normalized[] = new String[3];
  normalized = input.split(":");
  
  int mult = (normalized[0].contains("-")) ? -1 : 1;
  
  double degrees = Double.parseDouble(normalized[0]);
  
  double minToDeg = Double.parseDouble(normalized[1]) / 60.0;
  double secToDeg = Double.parseDouble(normalized[2]) / 3600.0;
  
  degrees += (minToDeg + secToDeg);
  
  degrees *= mult;
  
  double rads = ((Math.PI)/180.0) * degrees;
  
  return (float)rads;
}

public void parseCoreEvents(String[] data) {
  int set_len = data.length;
  events = new SNova[data.length];
  
  mean_v_band = 0.0;
  mean_i_band = 0.0;
  mean_r_band = 0.0;
  
  v_stats = new double[2];
  i_stats = new double[2];
  r_stats = new double[2];
  
  String partitionedLine[] = new String[numInitProperties];
  
  
  double max_v = Double.MIN_VALUE;
  double min_v = Double.MAX_VALUE;
  
  double max_i = Double.MIN_VALUE;
  double min_i = Double.MAX_VALUE;
  
  double max_r = Double.MIN_VALUE;
  double min_r = Double.MAX_VALUE;
  

  for (int i = 0 ; i < set_len; i++) {
    data[i] = data[i].replaceAll("\\s+", " ");
    partitionedLine = data[i].split(" ");
      
    String formattedRA = partitionedLine[1];
    String formattedDec = partitionedLine[2];
    float eventRA = getRA(formattedRA);
    float eventDec = getDec(formattedDec);
    String name = partitionedLine[0];
    String type = partitionedLine[3];
    String hostGalaxy = partitionedLine[4];
    
    float julianDate = Float.parseFloat(partitionedLine[5]);
    float v_band_mag = Float.parseFloat(partitionedLine[6]);
    float r_band_mag = Float.parseFloat(partitionedLine[10]);
    float i_band_mag = Float.parseFloat(partitionedLine[7]);
    
    
    if (max_v < v_band_mag) {
      max_v = v_band_mag;  
    }
    
    if (min_v > v_band_mag && v_band_mag != 0) {
      min_v = v_band_mag;  
    }
    
    if (max_i < i_band_mag) {
      max_i = i_band_mag;  
    }
    
    if (min_i > i_band_mag && i_band_mag != 0) {
      min_i = i_band_mag;  
    }
    
    if (max_r < r_band_mag) {
      max_r = r_band_mag;  
    }
    
    if (min_r > r_band_mag && r_band_mag != 0) {
      min_r = r_band_mag;
    }
    
    mean_r_band += r_band_mag;
    mean_v_band += v_band_mag;
    mean_i_band += i_band_mag;
    
    SNova event = new SNova(name, type, hostGalaxy, formattedRA, formattedDec, eventRA, eventDec);
    event.setCartesian(eventRA, eventDec, celestial_radius);
    event.setDateMagnitudes(julianDate, v_band_mag, r_band_mag, i_band_mag);
    //event.setGregorianDate(julianDate);
    events[i] = event;  
  }
 
  
  v_stats[0] = min_v;
  v_stats[1] = max_v;
  
  i_stats[0] = min_i;
  i_stats[1] = max_i;
  
  r_stats[0] = min_r;
  r_stats[1] = max_r;
  
  mean_r_band /= set_len;
  mean_v_band /= set_len;
  mean_i_band /= set_len;
  
  for (int i = 0; i < events.length; i ++) {
    events[i].computeColor();  
    
  }
}

public void parseBinaryEvents(String[] data) {
  int set_len = data.length;
  binary_events = new SNova[data.length];
  
  String partitionedLine[] = new String[numInitProperties];
  
  
  for (int i = 0 ; i < set_len; i++) {
    data[i] = data[i].replaceAll("\\s+", " ");
    partitionedLine = data[i].split(" ");
    
    String right_asc = partitionedLine[2] + ":" + partitionedLine[3] + ":" + 0;
    String dec = partitionedLine[4] + ":" + partitionedLine[5] + ":" + 0;
    
    float eventRA = getRA(right_asc);
    float eventDec = getDec(dec);
    String name = partitionedLine[0] + partitionedLine[1];
    String type = "Ia";
    String hostGalaxy = "N/A";
    float julianDate = Float.parseFloat(partitionedLine[6]);
    float v_band_mag = Float.parseFloat(partitionedLine[7]);

    SNova event = new SNova(name, type, hostGalaxy , right_asc, dec, eventRA, eventDec);
    event.setCartesian(eventRA, eventDec, celestial_radius);
    event.setDateMagnitudes(julianDate, v_band_mag, 0.0, 0.0);
    //event.setGregorianDate(julianDate);
    
    binary_events[i] = event; 
  }
}




void setup() {
  size(1000, 1000, P3D);
  //frame.addMouseWheelListener(new MouseWheelInput()); 
  data = loadFont("BodoniMT-14.vlw");
  textFont(data, 12);
  
  String sNovaInput[] = loadStrings("SNinfo.txt");
  String binNovaInput[] = loadStrings("binarySN.txt");
  
  parseCoreEvents(sNovaInput);
  parseBinaryEvents(binNovaInput);
  
  texmap = loadImage("mwpan2_Merc_2000x1200.jpg");    
  initializeSphere(sDetail);
  
    
}
 
void keyPressed() {
  switch(key) {
    case 'W':
    case 'w':
      delta_lon = 0.025;
      break;
    case'S':
    case 's':
      delta_lon = -0.025;
      break;
    case'A':
    case 'a':
      delta_lat = 0.025;
      break;
    case 'D':
    case 'd':
      delta_lat= -0.025;
      break;
    case ' ':
      lat = 0;
      lon = 0;
      break;
    default:
      break;
  }
 
}

void keyReleased() {
  switch(key) {
    case 'W':
    case 'w':
      delta_lon = 0.0;
      break;
    case'S':
    case 's':
      delta_lon = 0.0;
      break;
    case'A':
    case 'a':
      delta_lat = 0.0;
      break;
    case 'D':
    case 'd':
      delta_lat= 0.0;
      break;
    default:
      break;
  }
}



float zoom(float scale_val) {
  scale_val += mouseZoom / 50;
  scale_val = constrain(scale_val, 1.0, 2.0);
  
  return scale_val;
}

void draw() {
  background(0);
  stroke(255);
  renderGlobe();
  
  
}

void mouseWheel(MouseEvent event) {
  mouseZoom =  -event.getAmount();
  
}

/* These rendering functions are taken from the texturing example provided by Mike 'Flux' Chang, http://processing.org/examples/texturesphere.html */
void renderGlobe() {
  pushMatrix();
  translate(width * 0.5, height * 0.5, pushBack);
  pushMatrix();
  noFill();
  stroke(255,200);
  strokeWeight(2);
  smooth();
  popMatrix();
  lights();    
  pushMatrix();
  scale_val = zoom(scale_val);
  scale(scale_val);
  rotateX( radians(-rotationX) );  
  rotateY( radians(270 - rotationY) );
  fill(200);
  noStroke();
  textureMode(IMAGE);  
  texturedSphere(globeRadius, texmap);
  
   PVector event_vect;
  
  //End-point of edge
  PVector end_vect = new PVector();
  
  for (int ind = 0; ind < events.length; ind ++) {
    stroke(events[ind].rgb_colors[0], events[ind].rgb_colors[1], events[ind].rgb_colors[2]);
    event_vect = events[ind].vect;
    
    end_vect = event_vect.mult(event_vect, ((events[ind].v_band_mag) / (float) Math.log(events[ind].v_band_mag)) / 2, end_vect);
    
    line(event_vect.x, event_vect.y, event_vect.z, end_vect.x, end_vect.y, end_vect.z);
    events[ind].showMouseOver(end_vect);  
    
  }
  
  
  
  for (int ind = 0; ind < binary_events.length; ind ++) {
    stroke(255, 140, 0);
    event_vect = binary_events[ind].vect;
    
    end_vect = event_vect.mult(event_vect, ((binary_events[ind].v_band_mag) / (float) Math.log(binary_events[ind].v_band_mag)) / 2, end_vect);
    
    line(event_vect.x, event_vect.y, event_vect.z, end_vect.x, end_vect.y, end_vect.z);
    binary_events[ind].showMouseOver(end_vect);
    
  }
  
  popMatrix();  
  popMatrix();
  rotationX += velocityX;
  rotationY += velocityY;
  velocityX *= 0.95;
  velocityY *= 0.95;
  
  // Implements mouse control (interaction will be inverse when sphere is  upside down)
  if(mousePressed){
    velocityX += (mouseY-pmouseY) * 0.01;
    velocityY -= (mouseX-pmouseX) * 0.01;
  }
}

void initializeSphere(int res)
{
  sinLUT = new float[SINCOS_LENGTH];
  cosLUT = new float[SINCOS_LENGTH];

  for (int i = 0; i < SINCOS_LENGTH; i++) {
    sinLUT[i] = (float) Math.sin(i * DEG_TO_RAD * SINCOS_PRECISION);
    cosLUT[i] = (float) Math.cos(i * DEG_TO_RAD * SINCOS_PRECISION);
  }

  float delta = (float)SINCOS_LENGTH/res;
  float[] cx = new float[res];
  float[] cz = new float[res];
  
  // Calc unit circle in XZ plane
  for (int i = 0; i < res; i++) {
    cx[i] = -cosLUT[(int) (i*delta) % SINCOS_LENGTH];
    cz[i] = sinLUT[(int) (i*delta) % SINCOS_LENGTH];
  }
  
  // Computing vertexlist vertexlist starts at south pole
  int vertCount = res * (res-1) + 2;
  int currVert = 0;
  
  // Re-init arrays to store vertices
  sphereX = new float[vertCount];
  sphereY = new float[vertCount];
  sphereZ = new float[vertCount];
  float angle_step = (SINCOS_LENGTH*0.5f)/res;
  float angle = angle_step;
  
  // Step along Y axis
  for (int i = 1; i < res; i++) {
    float curradius = sinLUT[(int) angle % SINCOS_LENGTH];
    float currY = -cosLUT[(int) angle % SINCOS_LENGTH];
    for (int j = 0; j < res; j++) {
      sphereX[currVert] = cx[j] * curradius;
      sphereY[currVert] = currY;
      sphereZ[currVert++] = cz[j] * curradius;
    }
    angle += angle_step;
  }
  sDetail = res;
}

void texturedSphere(float r, PImage t) {
  int v1,v11,v2;
  r = (r + 240 ) * 0.33;
  beginShape(TRIANGLE_STRIP);
  texture(t);
  float iu=(float)(t.width-1)/(sDetail);
  float iv=(float)(t.height-1)/(sDetail);
  float u=0,v=iv;
  for (int i = 0; i < sDetail; i++) {
    vertex(0, -r, 0,u,0);
    vertex(sphereX[i]*r, sphereY[i]*r, sphereZ[i]*r, u, v);
    u+=iu;
  }
  vertex(0, -r, 0,u,0);
  vertex(sphereX[0]*r, sphereY[0]*r, sphereZ[0]*r, u, v);
  endShape();   
  
  // Middle rings
  int voff = 0;
  for(int i = 2; i < sDetail; i++) {
    v1=v11=voff;
    voff += sDetail;
    v2=voff;
    u=0;
    beginShape(TRIANGLE_STRIP);
    texture(t);
    for (int j = 0; j < sDetail; j++) {
      vertex(sphereX[v1]*r, sphereY[v1]*r, sphereZ[v1++]*r, u, v);
      vertex(sphereX[v2]*r, sphereY[v2]*r, sphereZ[v2++]*r, u, v+iv);
      u+=iu;
    }
  
    // Close each ring
    v1=v11;
    v2=voff;
    vertex(sphereX[v1]*r, sphereY[v1]*r, sphereZ[v1]*r, u, v);
    vertex(sphereX[v2]*r, sphereY[v2]*r, sphereZ[v2]*r, u, v+iv);
    endShape();
    v+=iv;
  }
  u=0;
  
  // Add the northern cap
  beginShape(TRIANGLE_STRIP);
  texture(t);
  for (int i = 0; i < sDetail; i++) {
    v2 = voff + i;
    vertex(sphereX[v2]*r, sphereY[v2]*r, sphereZ[v2]*r, u, v);
    vertex(0, r, 0,u,v+iv);    
    u+=iu;
  }
  vertex(sphereX[voff]*r, sphereY[voff]*r, sphereZ[voff]*r, u, v);
  endShape();
  
}