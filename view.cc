#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include <GL/glut.h>
#include <map>
#include <set>
#include <cmath>
#include "view.h"


GLdouble bodyWidth = 1.0;

GLfloat angle = -150;   /* in degrees */
GLfloat xloc = 0, yloc = 0, zloc = 0;
int moving, beginVic=0;
int newModel = 1;


/* ARGSUSED3 */
void mouse(int button, int state, int x, int y)
{
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    moving = 1;
    beginVic = x;
  }
  if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
    moving = 0;
  }
}

/* ARGSUSED1 */
void motion(int x, int y)
{
  if (moving) {
    angle = angle + (x - beginVic);
    beginVic = x;
    newModel = 1;
    glutPostRedisplay();
  }
}

void
tablet(int x, int y)
{
  xloc = ((GLfloat) x) / 500 - 4;
  yloc = ((GLfloat) y) / 1000 - 2;
  newModel = 1;
  glutPostRedisplay();
}

int xt = 1, yt = 1, zt = 1, xr = 1;

void 
translate(int x, int y, int z)
{
  GLfloat newz;

  if (xt)
    xloc += ((GLfloat) x) / 100;
  if (yt)
    yloc += ((GLfloat) y) / 100;
  if (zt) {
    newz = zloc - ((GLfloat) z) / 100;
    if (newz > -60.0 && newz < 13.0)
      zloc = newz;
  }
  newModel = 1;
  glutPostRedisplay();
}

/* ARGSUSED1 */
void
rotate(int x, int y, int z)
{
  if (xr) {
    angle += x / 2.0;
    newModel = 1;
    glutPostRedisplay();
  }
}

void
button(int button, int state)
{
  if (state == GLUT_DOWN) {
    switch (button) {
    case 1:
      xt = yt = zt = xr = 1;
      break;
    case 5:
      xt = 1;
      yt = zt = xr = 0;
      break;
    case 6:
      yt = 1;
      xt = zt = xr = 0;
      break;
    case 7:
      zt = 1;
      xt = yt = xr = 0;
      break;
    case 8:
      xr = 1;
      xt = yt = zt = 0;
      break;
    case 9:
      xloc = yloc = zloc = 0;
      newModel = 1;
      glutPostRedisplay();
      break;
    }
  }
}




int nRows = 480;
int nCols = 480; 

GLfloat light_ambient[] = {0.5, 0.5, 0.5, 1.0};  /* Red diffuse light. */
GLfloat light_diffuse[] = {0.8, 0.8, 0.8, 1.0};  /* Red diffuse light. */
GLfloat light_specular[] = {0.8, 0.8, 0.8, 1.0};  /* Red diffuse light. */
GLfloat light_position[] = {0.0, 0.0, 1.0, 0.0};  /* Infinite light location. */

static float modelAmb[4] = {0.2, 0.2, 0.2, 1.0};
static float matAmb[4] = {0.2, 0.2, 0.2, 1.0};
static float matDiff[4] = {0.8, 0.8, 0.8, 1.0};
static float matSpec[4] = {0.4, 0.4, 0.4, 1.0};
static float matEmission[4] = {0.0, 0.0, 0.0, 1.0};

static float modelAmb2[4] = {0.5, 0.5, 0.5, 1.0};
static float matAmb2[4] = {0.5, 0.5, 0.5, 1.0};
static float matDiff2[4] = {0.8, 0., 0., 1.0};
static float matSpec2[4] = {0.4, 0., 0., 1.0};
static float matEmission2[4] = {0.0, 0.0, 0.0, 1.0};



TriangleMesh trig;
GLUquadricObj *qobj;



/*
bool contain(Edge & e, map < pair <int, int> , Edge > & list) 
{

	pair <int, int> key;

	key.first = e.v1;
	key.second = e.v2;

	if (list.find(key) == list.end()) return false;
	else return true;
}
*/

bool contain(Edge & e, vector < Edge > & list) 
{
	bool inlist = false;

	for (int i = 0; i < list.size(); i++) 
	{
		if ((list[i]._v1 == e._v1 && list[i]._v2 == e._v2) ||
		    (list[i]._v2 == e._v1 && list[i]._v1 == e._v2)) 
		{
			return true;
		}	
	}

	return false;
}

int edgeID(Edge & e, vector < Edge > & list) 
{
	bool inlist = false;

	for (int i = 0; i < list.size(); i++) 
	{
		if ((list[i]._v1 == e._v1 && list[i]._v2 == e._v2) ||
		    (list[i]._v2 == e._v1 && list[i]._v1 == e._v2)) 
		{
			return i;
		}	
	}

	return -1;
}




/*
int edgeID(Edge & e, map < pair <int, int> , Edge > & list) 
{
	pair <int, int> key;

	key.first = e.v1;
	key.second = e.v2;

	if (list.find(key) == list.end()) return -1;
	else return list[key].id();   
}
*/



int find(Edge & e, vector <Edge> list) 
{
	for (int i = 0; i < list.size(); i++) {
		if (list[i] == e) return i;
	}

	return -1;
}


void TriangleMesh::loadFile(char * filename)
{
	ifstream f(filename);

	if (f == NULL) {
		cerr << "failed reading polygon data file " << filename << endl;
		exit(1);
	}

	char buf[1024];
	char header[100];
	float x,y,z;
	float xmax,ymax,zmax,xmin,ymin,zmin;
	int v1, v2, v3, n1, n2, n3;

	xmax =-10000; ymax =-10000; zmax =-10000;
	xmin =10000; ymin =10000; zmin =10000;

	while (!f.eof()) {
		    f.getline(buf, sizeof(buf));
		    sscanf(buf, "%s", header);  

		    if (strcmp(header, "v") == 0) {
				sscanf(buf, "%s %f %f %f", header, &x, &y, &z);  
				_v.push_back(Vector3f(x,y,z));

				_vn.push_back(Vector3f(0.f,0.f,1.f));

				Node node;

				node._id = _v.size()-1; 

				_node.push_back(node);
			

				if (x > xmax) xmax = x;
				if (y > ymax) ymax = y;
				if (z > zmax) zmax = z;

				if (x < xmin) xmin = x;
				if (y < ymin) ymin = y;
				if (z < zmin) zmin = z;
		    }
		    else if (strcmp(header, "vn") == 0) {
		//	sscanf(buf, "%s %f %f %f", header, &x, &y, &z);  
		//	_vn.push_back(Vector3f(x,y,z));
		    }
		    else if (strcmp(header, "f") == 0) 
		    {
		//	sscanf(buf, "%s %d//%d %d//%d %d//%d", header, &v1, &n1,
		//		&v2, &n2, &v3, &n3);
			

			sscanf(buf, "%s %d %d %d", header, &v1, &v2, &v3);


			Triangle trig(v1-1, v2-1, v3-1, v1-1, v2-1, v3-1);
			trig._id = _trig.size(); 
			_trig.push_back(trig);

			Edge e1(v1-1, v2-1);
			Edge e2(v2-1, v3-1);
			Edge e3(v3-1, v1-1);

			/*
			pair <int, int> id10(v1-1,v2-1),id11(v2-1,v1-1),
			     		id20(v2-1,v3-1),id21(v3-1,v2-1),
					id30(v3-1,v1-1),id31(v1-1,v3-1); 
			*/

			int id1,id2,id3;

			if ((id1 = edgeID(e1, _edge)) < 0) 
			{
//				_edge[id10] = e1; _edge[id11] = e1;

//				_edge[id10].setId(_edge.size()/2);
//				_edge[id11].setId(_edge.size()/2);

//				_edge[id10].add_triangle(&trig);
//				_edge[id11].add_triangle(&trig);

				id1 = _edge.size();
				_edge.push_back(e1);
				_edge[_edge.size()-1] = e1;

				_node[v1-1].edges_to.push_back(v2-1);
				_node[v2-1].edges_to.push_back(v1-1);


				_node[v1-1].edges_cost.push_back(_v[v1-1].distance(_v[v2-1]));
				_node[v2-1].edges_cost.push_back(_v[v1-1].distance(_v[v2-1]));
			}

			if ((id2 = edgeID(e2, _edge)) < 0) 
			{
				/*
				_edge[id20] = e2; _edge[id21] = e2;
				
				_edge[id20].setId(_edge.size()/2);
				_edge[id21].setId(_edge.size()/2);

				_edge[id20].add_triangle(&trig);
				_edge[id21].add_triangle(&trig);
				*/

				id2 = _edge.size();
				e2.setId(id2);
				e2.add_triangle(trig._id);
				_edge.push_back(e2);
				_edge[_edge.size()-1] = e2;

				_node[v2-1].edges_to.push_back(v3-1);
				_node[v3-1].edges_to.push_back(v2-1);


				_node[v2-1].edges_cost.push_back(_v[v2-1].distance(_v[v3-1]));
				_node[v3-1].edges_cost.push_back(_v[v3-1].distance(_v[v2-1]));
			}

			if ((id3 = edgeID(e3, _edge)) < 0) 
			{
				/*
				_edge[id30] = e3; _edge[id31] = e3;

				_edge[id30].setId(_edge.size()/2);
				_edge[id31].setId(_edge.size()/2);

				_edge[id30].add_triangle(&trig);
				_edge[id31].add_triangle(&trig);
				*/

				id3 = _edge.size();
				e3.setId(id3);
				e3.add_triangle(trig._id);
				_edge.push_back(e3);

				_node[v3-1].edges_to.push_back(v1-1);
				_node[v1-1].edges_to.push_back(v3-1);


				_node[v3-1].edges_cost.push_back(_v[v3-1].distance(_v[v1-1]));
				_node[v1-1].edges_cost.push_back(_v[v1-1].distance(_v[v3-1]));
			}

			_edge[id1].add_triangle(trig._id);
			_edge[id2].add_triangle(trig._id);
			_edge[id3].add_triangle(trig._id);


//			_trig[_trig.size()-1].setEdge(_edge[id10].id(), _edge[id20].id(), _edge[id30].id());
			_trig[_trig.size()-1].setEdge(id1,id2,id3); 
//			_trig[_trig.size()-1].setEdge(&_edge[id1], &_edge[id2], &_edge[id3]);

			//cout << " set Edge "<< "ids " << id1 << ' '<< id2 << ' '<<id3<<'-' << _edge[id1].id() << ' ' <<  _edge[id2].id() << ' ' <<  _edge[id3].id() << endl;
			/*
			int tmpid = _trig.size()-1 ;
			cout << " trig " << _trig.size()-1 << ' '; 
			cout << _trig[tmpid].edge(0) << ' ' << _trig[tmpid].edge(1) << ' ' 
			     << _trig[tmpid].edge(2) << endl; 
			*/

		    }
 	}

	vector < vector < int > > facelist (_v.size());
	vector < Vector3f > facenorm (_trig.size());

	for (int i = 0; i < _edge.size(); i++) {
		//cout << " edge " << i << " trig list " << _edge[i].getTrigList().size()<< endl;  //XXXXXXXXXXXX
	}

	for (int i = 0; i < _trig.size(); i++) 
	{
		/*
		cout << " trig " << i << ' '; 
			cout << _trig[i].edge(0) << ' ' << _trig[i].edge(1) << ' ' 
			     << _trig[i].edge(2) << endl; 
		*/


		Vector3f tmpv = (_v[_trig[i]._vertex[2]] - _v[_trig[i]._vertex[0]]) % 
				(_v[_trig[i]._vertex[1]] - _v[_trig[i]._vertex[0]]) ;

		tmpv.normalize();
		facenorm[i] = tmpv;

		facelist[_trig[i]._vertex[0]].push_back(i);
		facelist[_trig[i]._vertex[1]].push_back(i);
		facelist[_trig[i]._vertex[2]].push_back(i);
	}


	for (int i = 0; i < _v.size(); i++)  
	{
		Vector3f N(0.f,0.f,0.f); 

		float rate1, rate2;

		if (_v[i][1] > 0.5) 
		{
		       rate1 = 1.f ; rate2 = 0.f;	
		}
		else if (_v[i][1] < -0.5) 
		{
		       rate1 = 0.f ; rate2 = 1.f;	
		}
		else 
		{
			rate1 = _v[i][1] + 0.5f; rate2 = 1.f - rate1; 
		}

	//	cout << " v " << i << " 1:" << rate1 << " 2:" << rate2 << endl ;

		for (int j = 0; j < facelist[i].size(); j++) 
		{
			N += facenorm[facelist[i][j]]; 
	//		cout << " f " << facelist[i][j] << ' ' ;
		}
	//	cout << endl;

		N /= (float)facelist[i].size();

		_vn[i] = N;
	}


	_xmin = xmin; _ymin = ymin; _zmin = zmin;
	_xmax = xmax; _ymax = ymax; _zmax = zmax;

	f.close();

};


void recalcModelView(void)
{
	glPopMatrix();
	glPushMatrix();
	glTranslatef(xloc, yloc, zloc);
	glRotatef(angle, 0.0, 1.0, 0.0);
	glTranslatef(0, 0, .0);
	newModel = 0;
}

void myDisplay()
{


	if (newModel)
		recalcModelView();


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear OpenGL Window
	int trignum = trig.trigNum();
	Vector3f v1,v2,v3,n1,n2,n3;


	for (int i = 0 ; i < trignum; i++)  
	{
		/*** do the rasterization of the triangles here using glRecti ***/
//		glColor3f(1,0,0);/* you can set the color by doing the computation of illumination*/
//		glRecti(0,0,1,1);/* filling in the pixel */

		float m1,m2,m3,min,max;
		trig.getTriangleVertices(i,v1,v2,v3);
		trig.getTriangleNormals(i,n1,n2,n3);
		trig.getMorseValue(i, m1, m2, m3);

		m1 = m2 = m3 = trig.color(i);

		GLfloat skinColor[] = {0.1, 1., 0.1, 1.0};

		if (max >= 0) {
			glBegin(GL_TRIANGLES);

				skinColor[1] = m1; skinColor[0] = 1-m1;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n1[0],-n1[1],-n1[2]);
				glVertex3f(v1[0],v1[1],v1[2]);

				skinColor[1] = m2; skinColor[0] = 1-m2;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n2[0],-n2[1],-n2[2]);
				glVertex3f(v2[0],v2[1],v2[2]);

				skinColor[1] = m3; skinColor[0] = 1-m3;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n3[0],-n3[1],-n3[2]);
				glVertex3f(v3[0],v3[1],v3[2]);

				skinColor[1] = m1; skinColor[0] = 1-m1;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n1[0],-n1[1],-n1[2]);
				glVertex3f(v1[0],v1[1],v1[2]);

			glEnd();
		}
	}

	 glutSwapBuffers();
}


//--------------------Write by Pengfei Gao------------------

#define PI 3.1415926
#define MAX_CHILDCOUNT	10
#define JOINT_NUM 22
#define JOINT_LINE 5
#define INTERPOLATION 18



typedef struct _Joint
{
	int id;
	int childCount;
	struct _Joint *child[MAX_CHILDCOUNT], *parent; //previous joint id
}Joint;

typedef struct _Frame
{
	vector<Vector3f> *translations;
	vector<Vector3f> *rotations;
}Frame;

vector<Edge> bone_list;
Joint *root;
vector<Vector3f> joint_list;
vector<Vector3f> trig_vback;

//The tree struct is used to store the structure of the skeleton
Joint *jointFindTree(Joint *root, int id)
{
	if(!root)
		return NULL;
	if(root->id == id)
		return root;
	Joint *res = (Joint *)malloc(sizeof(Joint));
	for (int i=0;i<root->childCount;i++)
		res = jointFindTree(root->child[i],id);
	return res;
}

Joint *jointAddChild(Joint *root, float id, float pid)
{
	Joint *t, *parent;
	
	if(!root){ 
		if(!(root = (Joint *)malloc(sizeof(Joint))))
			return NULL;
		root->parent = NULL;
		parent = root;
	}
	/* allocate the child */
	else if (root->childCount < MAX_CHILDCOUNT)
	{
		if((int)pid==-1)
			parent = root;
		else
			parent = jointFindTree(root,pid);
		if(parent == NULL)
			return NULL;
	
		if(!(t = (Joint *)malloc(sizeof(Joint))))
			return NULL; //error
		t->parent = parent; //set it's previous joint
		parent->child[parent->childCount++] = t;
		parent = t; //change the root
	}
	else //can't add a child
		return NULL;
	parent->id = id;
	parent->childCount = 0;
	
	for(int i=0;i<MAX_CHILDCOUNT;i++)
		parent->child[i] = NULL;
		
	return root;
}



Joint *readSkeleton(char* filename, vector<Vector3f> &joint_list)
{
	ifstream f(filename);
	if (f == NULL) {
		cerr << "failed reading skeleton data file " << filename << endl;
		exit(1);
	}
	
	char buf[1024];
	int i = 0;
	float data[22][5];
	while (!f.eof())
	{
		f.getline(buf, sizeof(buf));
		sscanf(buf, "%f %f %f %f %f", &data[i][0],&data[i][1],&data[i][2],&data[i][3],&data[i][4]);	
		i++;
	}
	f.close();

	for(int i=21;i>0;i--){
		int pid = data[i][4];
		data[i][1] = data[i][1] - data[pid][1];
		data[i][2] = data[i][2] - data[pid][2];
		data[i][3] = data[i][3] - data[pid][3];
	}

	Joint *root, *tmp;
	if(!(root = jointAddChild(NULL,data[0][0],data[0][4])))
	{
		fprintf(stderr, "Error! Can't create a root!\n");
		exit(EXIT_FAILURE);
	}
	Vector3f v(data[0][1],data[0][2],data[0][3]);
	joint_list.push_back(v);
	for(i=1;i<JOINT_NUM;i++)
	{
		root = jointAddChild(root,data[i][0],data[i][4]);
		Vector3f v(data[i][1],data[i][2],data[i][3]);
		joint_list.push_back(v);
	}
	
	return root;
	
}



vector<vector<float> > weights_list;
void readWeights(char* filename,vector<vector<float> > &weightslist)
{
	ifstream f(filename);
	
	string line;
	float data[JOINT_NUM][JOINT_LINE];
	while (getline(f, line))
	{
		istringstream ss(line);

		istream_iterator<float> beginVic(ss), end;
		vector<float> list(beginVic,end);
		weights_list.push_back(list);
	}
	f.close();
}

//read keyframes from file, one keyframe is 22 lines and each line contain translation coordinates and rotation angles
vector<Frame> keyframes;
void readKeyframes(char *filename, vector<Frame> &keyframes)
{
	ifstream f(filename);
	char buf[1024];
	float tx,ty,tz,rx,ry,rz;

	int i=0;
	Frame frame;
	frame.translations = new vector<Vector3f>();
	frame.rotations = new vector<Vector3f>();
	
	while (!f.eof()) {
		
		f.getline(buf, sizeof(buf));
		sscanf(buf, "%f %f %f %f %f %f %f", &tx, &ty, &tz, &rx, &ry, &rz); 
			
		frame.translations->push_back(Vector3f(tx,ty,tz));
		frame.rotations->push_back(Vector3f(rx,ry,rz));

		tx = ty = tz = rx = ry = rz = 0;
		i++;
		if(i==22){
			vector<Vector3f> *t = new vector<Vector3f>(*(frame.translations));
			vector<Vector3f> *r = new vector<Vector3f>(*(frame.rotations));
			
			Frame tmp;
			tmp.translations = t;
			tmp.rotations = r;

			keyframes.push_back(tmp);
			frame.translations->clear();
			frame.rotations->clear();
			i = 0;	
		}
			
	}
	f.close();
}


//interpolate 18 frames between two keyframes
vector<Frame> allframes; 
void interpolateFrames(vector<Frame> keyframes, vector<Frame> &allframes)
{
	int k = 0;
	Vector3f t,r;

	allframes.clear();
	for(int i=1;i<keyframes.size();i++){
		for(int k=0; k<INTERPOLATION;k++){
			Frame frame;
			frame.translations = new vector<Vector3f>();
			frame.rotations = new vector<Vector3f>();
			for(int j=0;j<JOINT_NUM;j++){
				t = keyframes[i].translations->at(j);
				t -= keyframes[i-1].translations->at(j);
				t /= INTERPOLATION;
				t *= k;
				t += keyframes[i-1].translations->at(j);
				r = keyframes[i].rotations->at(j);
				r -= keyframes[i-1].rotations->at(j);
				r /= INTERPOLATION;
				r *= k;
				r += keyframes[i-1].rotations->at(j);
				frame.translations->push_back(t);
				frame.rotations->push_back(r);
			} 
			allframes.push_back(frame);
		}
	}
}

//find all the bones and store into the bone_list
void skeletonEdge(Joint *root, vector<Edge> &bone_list)
{
	
	if(!root)
		return ;
	
	Joint *child;
	for(int i=0;i<root->childCount;i++)
	{
		child = root->child[i];
		Edge e(root->id, child->id);
		bone_list.push_back(e);
		skeletonEdge(child,bone_list);
	}
}


void translateVertex(Vector3f &v,Vector3f t)
{ 
	v[0] = v[0] + t[0];
	v[1] = v[1] + t[1];
	v[2] = v[2] + t[2];
}

void ratateXVertex(Vector3f &v, float theta)
{
	float a = v[1];
	float b = v[2];
	v[1] = a * cos(theta) - b * sin(theta);
	v[2] = a * sin(theta) + b * cos(theta);
}

void rotateYVertex(Vector3f &v, float theta)
{
	float a = v[0];
	float c = v[2];
	v[0] = a * cos(theta) + c * sin(theta);
	v[2] = - a * sin(theta) + c * cos(theta);
}

void rotateZVertex(Vector3f &v, float theta)
{
	float a = v[0];
	float b = v[1];
	v[0] = a * cos(theta) - b * sin(theta);
	v[1] = a * sin(theta) + b * cos(theta);
}

void rotation(Vector3f &v, Vector3f r)
{
	rotateZVertex(v, r[2] * PI / 180);
	rotateYVertex(v, r[1] * PI / 180);
	ratateXVertex(v, r[0] * PI /180);
}

void transrotation(Vector3f &v, Vector3f t, Vector3f r)
{
	translateVertex(v, t);
	rotation(v, r);
}

Vector3f inverseMatrix(const Vector3f v)
{
	Vector3f inv = v;
	inv *= -1;
	return inv;
}


//convert the local coordinates to world coordinates
void calcGlobalCoor(Frame frame,vector<Vector3f> &v)
{
	//vector<Vector3f> v(JOINT_NUM,Vector3f(0,0,0));
	Vector3f xx(0,0,0);
	for(int i=0;i<v.size();i++){
		v[i] = xx;
	}
	//root
	transrotation(v[0],frame.translations->at(0),frame.rotations->at(0));
	int link[5][8] = {
			{7,6,5,0},//head and spine
			{4,3,2,1,0}, //left leg
			{21,20,19,18,0}, //right leg
			{12,11,10,9,8,6,5,0}, //right arm
			{17,16,15,14,13,6,5,0}//left arm
			};
	vector<int> processed;
	int flag = 0;
	for (int i=0;i<5;i++){
		for(int j=0;j<8;j++){
			if(link[i][j]==0) break;
			for(int t=0;t<processed.size();t++){
				if(link[i][j] == processed.at(t)){ 
					flag = 1;
					break;
				}
			}
			if(flag){
				flag = 0;
				break;
			}
			for(int k=j;k<8;k++){
				transrotation( v[ link[i][j] ], frame.translations->at(link[i][k]),frame.rotations->at(link[i][k]) );
				if(link[i][k]==0) break;
			}
			processed.push_back(link[i][j]);
		}
	}	

}




vector<Vector3f> joint_world_coor(JOINT_NUM,Vector3f(0,0,0));

//draw the animations
void animateDisplay()
{
	if (newModel)
		recalcModelView();
		
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	

	glColor3f(1,0,0);
	int trignum = trig.trigNum();
	Vector3f v1,v2,v3,n1,n2,n3;

	for (int i = 0 ; i < trignum; i++)  
	{
		float m1,m2,m3,min,max=0;
		trig.getTriangleVertices(i,v1,v2,v3);
		trig.getTriangleNormals(i,n1,n2,n3);
		trig.getMorseValue(i, m1, m2, m3);

		m1 = m2 = m3 = trig.color(i);

		GLfloat skinColor[] = {0.1, 1.0, 0.1, 1.0};

		if (max >= 0) {
			glBegin(GL_TRIANGLES);

				skinColor[1] = m1; skinColor[0] = 1-m1;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n1[0],-n1[1],-n1[2]);
				glVertex3f(v1[0],v1[1],v1[2]);

				skinColor[1] = m2; skinColor[0] = 1-m2;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n2[0],-n2[1],-n2[2]);
				glVertex3f(v2[0],v2[1],v2[2]);

				skinColor[1] = m3; skinColor[0] = 1-m3;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n3[0],-n3[1],-n3[2]);
				glVertex3f(v3[0],v3[1],v3[2]);

				skinColor[1] = m1; skinColor[0] = 1-m1;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor); 
				glNormal3f(-n1[0],-n1[1],-n1[2]);
				glVertex3f(v1[0],v1[1],v1[2]);

			glEnd();
		}
	}

 	glutSwapBuffers();

}

void reshape (int w, int h)
{
	glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective(60.0, (GLfloat) w/(GLfloat) h, 1.0, 20.0);
	glMatrixMode (GL_MODELVIEW);
}

//respond to the keyboard events 
#include <pthread.h>
bool update=true;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
void keyboard(unsigned char key, int x, int y)
{
	pthread_mutex_lock( &mutex1 );
	switch (key)
	{
		case 27: //Escape key
			exit(0);break;
		case 32: //space
			update=false; break;
		case 119: //w
			update=false; 
			keyframes.clear();
			readKeyframes(const_cast<char *>("walk.txt"),keyframes); 
			allframes.clear();
			interpolateFrames(keyframes, allframes);
			update=true;
			break;
		case 115: //s
			update=false; 
			keyframes.clear();
			readKeyframes(const_cast<char *>("skip.txt"),keyframes);  
			allframes.clear();
			interpolateFrames(keyframes, allframes);
			update=true;
			break;
		case 97: //a
			update=false; 
			keyframes.clear();
			readKeyframes(const_cast<char *>("attack.txt"),keyframes);  
			allframes.clear();
			interpolateFrames(keyframes, allframes);
			update=true;
			break;
		case 100: //d
			update=false; 
			keyframes.clear();
			readKeyframes(const_cast<char *>("defend.txt"),keyframes);  
			allframes.clear();
			interpolateFrames(keyframes, allframes);
			update=true;
			break;
	}
	pthread_mutex_unlock( &mutex1 );
}

int i = 1;
void animate(int value)
{	
	if(update){
		vector<Vector3f> last_coor(JOINT_NUM,Vector3f(0,0,0));
		calcGlobalCoor(allframes[0], last_coor);
		calcGlobalCoor(allframes[i],joint_world_coor);
		i++;
		if(i==allframes.size())
			i=1;

		Vector3f temp(0,0,0);
		for(int i=0;i<trig._v.size();i++){
			for(int j=0;j<21;j++){
				if(weights_list[i][j] != 0)	{
					Vector3f bone_last(0,0,0);
					Vector3f bone_current(0,0,0);
					bone_last += last_coor[bone_list[j].v1()];
					bone_last += last_coor[bone_list[j].v2()];
					bone_last /= 2;
					bone_current += joint_world_coor[bone_list[j].v1()];		
					bone_current += joint_world_coor[bone_list[j].v2()];
					bone_current /= 2;
					temp[0] += weights_list[i][j] * ( trig_vback[i][0] + inverseMatrix(bone_last)[0] + bone_current[0] );
					temp[1] += weights_list[i][j] * ( trig_vback[i][1] + inverseMatrix(bone_last)[1] + bone_current[1] );
					temp[2] += weights_list[i][j] * ( trig_vback[i][2] + inverseMatrix(bone_last)[2] + bone_current[2] );
				}	
			}
			trig._v[i] = temp;
			temp -= temp;
		}
	}

	glutPostRedisplay();
	glutTimerFunc(24,animate,0);
}

//--------------------end------------------



int main(int argc, char **argv)
{
	if (argc >  1)  {
		trig.loadFile(argv[1]);
	}
	else {
		cerr << argv[0] << " <filename> " << endl;
		exit(1);
	}
	
	trig_vback.clear();
	for(int i=0;i<trig._v.size();i++)
		trig_vback.push_back(trig._v[i]);

	root = readSkeleton(const_cast<char *>("skeleton2.out"),joint_list);
	skeletonEdge(root, bone_list);
	readWeights(const_cast<char *>("attachment2.out"),weights_list);

	readKeyframes(const_cast<char *>("motion.txt"),keyframes);
	interpolateFrames(keyframes, allframes);

	calcGlobalCoor(allframes[0], joint_world_coor);

	
	int width, height;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);

	glutInitWindowSize(nRows, nCols);
	glutCreateWindow("SimpleExample");

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);


	/* Use depth buffering for hidden surface elimination. */
	glEnable(GL_DEPTH_TEST);

	/* Setup the view of the cube. */
	glMatrixMode(GL_PROJECTION);
	gluPerspective( /* field of view in degree */ 40.0, /* aspect ratio */ 1., /* Z near */ 1.0, /* Z far */ 1000.0);

	glMatrixMode(GL_MODELVIEW);

	gluLookAt(0.0, 0.0, 7.0,  /* eye is at (0,0,5) */
		  0.0, 0.0, 0.0,      /* center is at (0,0,0) */
		  0.0, 1.0, 0.0);      /* up is in positive Y direction */
	glPushMatrix();       /* dummy push so we can pop on model recalc */

	glutDisplayFunc(animateDisplay);// Callback function
	glutReshapeFunc(reshape);
 
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard); //control keyboard
	glutMotionFunc(motion);
	glutTabletMotionFunc(tablet);
	glutSpaceballMotionFunc(translate);
	glutSpaceballRotateFunc(rotate);
	glutSpaceballButtonFunc(button);
	
	
	glutTimerFunc(24,animate,0);

	glutMainLoop();// Display everything and wait
	
	return 0;

}

