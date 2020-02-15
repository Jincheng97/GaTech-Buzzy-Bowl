/*
Author: Jincheng Zhu
Class: ECE6122
Last Modified Date: 12/3/2019
Description: This is the code of final project.
First, I let the UAVs fly directly towards the center of the sphere (accelerate, fly with a constant speed and decelerate).
Then, I let the UAVs fly towards the equator of the sphere.
Finally, I let the UAVs fly along the meridians of the sphere.
Since the UAVs are not actually flying along the circular curves on the surface of the sphere, it is possible that the UAVs 
somehow fly just beneath the surface of the sphere. Thus, I made some proper approximation of the radius. 
Also, since the UAVs are able to generate a single force vector with a total magnitude of 20 N in any direction, while the 
gravity of each UAV is 10 N, I can consider a situation where the UAVs are flying in a environment with no gravity, and the
force they can generate is 10 N.
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "iomanip"
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <chrono>
#include <thread>
#include "ECE_Bitmap.h"


#define PI 3.1415926

typedef struct Image Image;
GLuint texture[1];
BMP inBitmap;

float red = 255;	// This is the red color of the UAVs.
int colorChange = 0;	// This helps counting the time for changing the color.

void drawFootballField(GLuint texture[])
{
	// Draw football field.
	glPushMatrix();
	glColor3f(0.0, 1.0, 0.0);
	
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glBegin(GL_QUADS);
	glTexCoord2f(1.0, 1.0);
	glVertex3f(-30, -60, 0.0);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(-30, 60, 0.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(30, 60, 0.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(30, -60, 0.0);
	glEnd();
	glBindTexture(GL_TEXTURE_2D, 0);
	glPopMatrix();

	glPushMatrix();
	glColor3f(0.5, 0.5, 0.5);
	glTranslatef(0.0, 0.0, 50.0);
	glutWireSphere(10, 10, 10);
	glPopMatrix();

}

void initTexture()
{
	glClearColor(0.5, 0.5, 0.5, 0.0);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	inBitmap.read("AmFBfield.bmp");

	glGenTextures(1, texture);

	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, inBitmap.bmp_info_header.width, inBitmap.bmp_info_header.height, 0, GL_RGB, GL_UNSIGNED_BYTE, &inBitmap.data[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glEnable(GL_TEXTURE_2D);
}

// This function draws all UAVs and the football field.
void drawUAVs()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawFootballField(texture);

	// Draw UAVs.
	float* sendBuf;	// This is the send buffer of MPI cimmunication.
	float* recvBuf;	// This is the receive buffer of MPI cimmunication.
	// These don't really matter. Just to fill in the blank units in the send and receive buffer.
	float UAVLocX = 10, UAVLocY = 10, UAVLocZ = 10;
	float UAVSpeedX = 0, UAVSpeedY = 0, UAVSpeedZ = 0;

	// Gather the location and speed info of all UAVs.
	sendBuf = new float[6];
	sendBuf[0] = UAVLocX;
	sendBuf[1] = UAVLocY;
	sendBuf[2] = UAVLocZ;
	sendBuf[3] = UAVSpeedX;
	sendBuf[4] = UAVSpeedY;
	sendBuf[5] = UAVSpeedZ;
	recvBuf = new float[16 * 6];
	MPI_Allgather(sendBuf, 6, MPI_FLOAT, recvBuf, 6, MPI_FLOAT, MPI_COMM_WORLD);
	
	glColor3f(red / 255, 0.0, 0.0);
	for (int i = 1; i < 16; i++)
	{
		glPushMatrix();
		glTranslatef(recvBuf[i * 6], recvBuf[i * 6 + 1], recvBuf[i * 6 + 2]);
		glutSolidTeapot(1);
		glPopMatrix();
	}
	glutSwapBuffers();
}

// This function sets the view port and perspective.
void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);	// This sets the view port.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (GLfloat)w / (GLfloat)h, 0.1, 1000.0);	// This sets the prospective.
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(100, 0, 100, 0.0, 0.0, 20.0, 0.0, 0.0, 1.0);
	drawFootballField(texture);
	glutSwapBuffers();
}

// This function is the timer.
void timer(int id)
{
	// The color of all UAVs ranges from 128 ~ 255 every 100 msc (according to our professor in class).
	if (colorChange % 256 < 127)
		red--;
	else
		red++;
	colorChange++;
	glutPostRedisplay();
	glutTimerFunc(100, timer, 0);
}

// This initializes OpenGL.
void mainOpenGL(int argc, char **argv)
{
    int mode = GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH;
    glutInit(&argc, argv);
    glutInitDisplayMode(mode);
    glutInitWindowSize(400, 400);
    glutCreateWindow("Football field");

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	initTexture();
    glutDisplayFunc(drawUAVs);
	glutReshapeFunc(reshape);

	glutTimerFunc(100, timer, 0);
	glutMainLoop();
}

// This is the main function.
int main(int argc, char** argv)
{
	// Initiate MPI.
	int numTasks, rank;
	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS)
	{
		printf("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	float UAVLocX, UAVLocY, UAVLocZ;	// These are the coordinates of the UAV.
	float UAVSpeedX, UAVSpeedY, UAVSpeedZ;	// These are the speed of the UAV.
	float* sendBuf;	// This is the send buffer used in MPI communication.
	float* recvBuf;	// This is the receive buffer used in MPI communication.
	float clockTime = 0;	// This is the clock.
	float dist;		// This is the distance from each UAV to the center of the sphere.
	float speed;	// This is the orbiting speed of each UAV.
	float time0, time1, time2;	// This is the time required in each stage.
	int state = 1;	// This helps deciding whether a UAV has flied pass the equator.
	float theta = 0;	// This is the radius angle.
	float ratioXY = 0;	// This is the ratio between the first two coordinates of a UAV.
	float signX = 0;	// This is the sign of the first coordinate.
	float signY = 0;	// This is the sign of the second coordinate.
	int count = 0;	// This helps counting the time.
	int inside = 0;	// This helps deciding whether a UAV have flied into the sphere with a radius of 10.

	// Set the orbiting speed of all UAVs (ranging from 5 ~ 10).
	speed = (rank - 1) / 3 + 5;

	// Set the initial location of all UAVs.
	if (rank <= 5)
		UAVLocX = -25;
	else if (rank <= 10)
		UAVLocX = 0;
	else if (rank <= 15)
		UAVLocX = 25;
	if (rank % 5 == 1)
		UAVLocY = -50;
	else if (rank % 5 == 2)
		UAVLocY = -25;
	else if (rank % 5 == 3)
		UAVLocY = 0;
	else if (rank % 5 == 4)
		UAVLocY = 25;
	else if (rank % 5 == 0)
		UAVLocY = 50;
	UAVLocZ = 0;

	// Set the initial speed of all UAVs.
	UAVSpeedX = 0;
	UAVSpeedY = 0;
	UAVSpeedZ = 0;
	
	// Calculate the distance from each UAV to the center of sphere.
	dist = sqrt(pow(UAVLocX, 2) + pow(UAVLocY, 2) + pow(UAVLocZ - 50, 2));

	// Calculate the time required in each stage. 
	time0 = 2 / 10;		// The time of acceleration and deceleration before a UAV reaches the surface of the sphere.
	time1 = (dist - 10.7 - 10 * pow(time0, 2)) / 2;		/* The time of a UAV flying with constant speed towards the surface of the sphere. 
														Since the flying routes of UAVs are not real circular curves, I made some approximation of the radius (10 -> 10.7).*/
	time2 = speed / 10;		// The time required for a UAV to reach the equator.

	if (rank == 0)	// The main thread runs OpenGL.
		mainOpenGL(argc, argv);
	else  // Other threads calculate the speed and location of UAVs.
	{
		while (clockTime < 60)	// The UAVs will fly for 60 sec.
		{
			float dirX, dirY, dirZ;	// These are the direction vector from the UAV to the center of sphere.
			float dirX0 = 0, dirY0 = 0, dirZ0 = 0;	// These are the vector perpendicular to the direction vector.
			float dir, dir0;	// These help normalizing the vectors.
			
			if (count % 10 == 0)	// The main thread will redisplay all UAVs every 0.01*10 seconds (100 msc).
			{
				sendBuf = new float[6];
				sendBuf[0] = UAVLocX;
				sendBuf[1] = UAVLocY;
				sendBuf[2] = UAVLocZ;
				sendBuf[3] = UAVSpeedX;
				sendBuf[4] = UAVSpeedY;
				sendBuf[5] = UAVSpeedZ;
				recvBuf = new float[16 * 6];
				MPI_Allgather(sendBuf, 6, MPI_FLOAT, recvBuf, 6, MPI_FLOAT, MPI_COMM_WORLD);

				// Compute whether two UAVs have collided. If they collide, they exchange the velocity vector.
				for (int i = 1; i < 15; i++)
				{
					for (int j = i + 1; j < 16; j++)
					{
						float UAVDist;
						UAVDist = sqrt(pow(recvBuf[i * 6] - recvBuf[j * 6], 2)
							+ pow(recvBuf[i * 6 + 1] - recvBuf[j * 6 + 1], 2)
							+ pow(recvBuf[i * 6 + 2] - recvBuf[j * 6 + 2], 2));
						if (UAVDist < 0.01)
						{
							float temp;
							for (int k = 0; k < 6; k++)
							{
								temp = recvBuf[i * 6 + k];
								recvBuf[i * 6 + k] = recvBuf[j * 6 + k];
								recvBuf[j * 6 + k] = temp;
							}
						}
					}
				}

				UAVLocX = recvBuf[rank * 6];
				UAVLocY = recvBuf[rank * 6 + 1];
				UAVLocZ = recvBuf[rank * 6 + 2];
				UAVSpeedX = recvBuf[rank * 6 + 3];
				UAVSpeedY = recvBuf[rank * 6 + 4];
				UAVSpeedZ = recvBuf[rank * 6 + 5];

				/* Check whether the UAVs have come within 10m of the center of the sphere.
				(The UAVs are not actually flying along the sphere, since they are taking the short cuts (circular archs).
				Thus, I set the thread hold as 9.3 instead of 10.) */
				for (int i = 1; i < 16; i++)
					if (sqrt(pow(recvBuf[i * 6], 2) + pow(recvBuf[i * 6 + 1], 2) + pow(recvBuf[i * 6 + 2] - 50, 2)) < 9.3)
						inside++;
				if (inside == 15)
					return 0;
			}

			// First the UAVs stay on the ground for 5 sec.
			if (count == 0)
				std::this_thread::sleep_for(std::chrono::seconds(5));
			// First let each UAV accelerates to 2 m/s and flies towards the sphere. 
			if (clockTime < time0)
			{
				// Calculate the flying direction of each UAV.
				dirX = 0 - UAVLocX;
				dirY = 0 - UAVLocY;
				dirZ = 50 - UAVLocZ;
				dir = sqrt(pow(dirX, 2) + pow(dirY, 2) + pow(dirZ, 2));
				dirX /= dir;
				dirY /= dir;
				dirZ /= dir;
				// Calculate the speed and location of each UAV.
				UAVLocX += UAVSpeedX * 0.01 + 1 / 2 * 10 * dirX * pow(0.01, 2);
				UAVLocY += UAVSpeedY * 0.01 + 1 / 2 * 10 * dirY * pow(0.01, 2);
				UAVLocZ += UAVSpeedZ * 0.01 + 1 / 2 * 10 * dirZ * pow(0.01, 2);
				UAVSpeedX += 10 * dirX * 0.01;
				UAVSpeedY += 10 * dirY * 0.01;
				UAVSpeedZ += 10 * dirZ * 0.01;
			}
			// Then let each UAV fly with a constant speed (2 m/s) towards the sphere.
			else if (clockTime < time0 + time1)
			{
				dirX = 0 - UAVLocX;
				dirY = 0 - UAVLocY;
				dirZ = 50 - UAVLocZ;
				dir = sqrt(pow(dirX, 2) + pow(dirY, 2) + pow(dirZ, 2));
				dirX /= dir;
				dirY /= dir;
				dirZ /= dir;
				UAVLocX += UAVSpeedX * 0.01;
				UAVLocY += UAVSpeedY * 0.01;
				UAVLocZ += UAVSpeedZ * 0.01;
				UAVSpeedX = 2 * dirX;
				UAVSpeedY = 2 * dirY;
				UAVSpeedZ = 2 * dirZ;
			}
			// Then let each UAV decelerate to 0 m/s and keeps flying towards the sphere.
			else if (clockTime < 2 * time0 + time1)
			{
				dirX = 0 - UAVLocX;
				dirY = 0 - UAVLocY;
				dirZ = 50 - UAVLocZ;
				dir = sqrt(pow(dirX, 2) + pow(dirY, 2) + pow(dirZ, 2));
				dirX /= dir;
				dirY /= dir;
				dirZ /= dir;
				UAVLocX += UAVSpeedX * 0.01 - 1 / 2 * 10 * dirX * pow(0.01, 2);
				UAVLocY += UAVSpeedY * 0.01 - 1 / 2 * 10 * dirY * pow(0.01, 2);
				UAVLocZ += UAVSpeedZ * 0.01 - 1 / 2 * 10 * dirZ * pow(0.01, 2);
				UAVSpeedX -= 10 * dirX * 0.01;
				UAVSpeedY -= 10 * dirY * 0.01;
				UAVSpeedZ -= 10 * dirZ * 0.01;
			}
			// Then let each UAV accelerate to the speed set before along the sphere surface.
			else if (clockTime < 2 * time0 + time1 + time2 - 0.03)
			{
				dirX = UAVLocX;
				dirY = UAVLocY;
				dirZ = UAVLocZ - 50;
				dir = sqrt(pow(dirX, 2) + pow(dirY, 2) + pow(dirZ, 2));
				dirX /= dir;
				dirY /= dir;
				dirZ /= dir; 
				if (dirX > 0)
					dirX0 = rand() % 600 / 100 + 1;
				else
					dirX0 = -(rand() % 600 / 100 + 1);
				if (dirY > 0)
					dirY0 = rand() % 600 / 100 + 1;
				else
					dirY0 = -(rand() % 600 / 100 + 1);
				dirZ0 = -(dirX * dirX0 + dirY * dirY0) / dirZ;
				dir0 = sqrt(pow(dirX0, 2) + pow(dirY0, 2) + pow(dirZ0, 2));
				dirX0 /= dir0;
				dirY0 /= dir0;
				dirZ0 /= dir0;

				UAVLocX += UAVSpeedX * 0.01 + 1 / 2 * 10 * dirX0 * pow(0.01, 2);
				UAVLocY += UAVSpeedY * 0.01 + 1 / 2 * 10 * dirY0 * pow(0.01, 2);
				UAVLocZ += UAVSpeedZ * 0.01 + 1 / 2 * 10 * dirZ0 * pow(0.01, 2);
				UAVSpeedX += 10 * dirX0 * 0.01;
				UAVSpeedY += 10 * dirY0 * 0.01;
				UAVSpeedZ += 10 * dirZ0 * 0.01;
			}
			// Finally let the UAVs fly towards the equator and fly along the meridians around the sphere.
			else
			{
				// First let the UAVs fly towards the equator.
				if (state == 1)
				{
					dirX = UAVLocX;
					dirY = UAVLocY;
					dirZ = UAVLocZ - 50;
					dir = sqrt(pow(dirX, 2) + pow(dirY, 2) + pow(dirZ, 2));
					dirX /= dir;
					dirY /= dir;
					dirZ /= dir;
					if (dirX > 0)
						dirX0 = rand() % 600 / 100 + 1;
					else
						dirX0 = -(rand() % 600 / 100 + 1);
					if (dirY > 0)
						dirY0 = rand() % 600 / 100 + 1;
					else
						dirY0 = -(rand() % 600 / 100 + 1);
					dirZ0 = -(dirX * dirX0 + dirY * dirY0) / dirZ;

					dir0 = sqrt(pow(dirX0, 2) + pow(dirY0, 2) + pow(dirZ0, 2));
					dirX0 /= dir0;
					dirY0 /= dir0;
					dirZ0 /= dir0;

					UAVLocX += UAVSpeedX * 0.01;
					UAVLocY += UAVSpeedY * 0.01;
					UAVLocZ += UAVSpeedZ * 0.01;
					UAVSpeedX = speed * dirX0;
					UAVSpeedY = speed * dirY0;
					UAVSpeedZ = speed * dirZ0;

					// Check if a UAV has reached the equator.
					if (fabs(dirZ0) > 0.9999)
					{
						signX = fabs(UAVLocX) / UAVLocX;
						signY = fabs(UAVLocY) / UAVLocY;
						ratioXY = fabs(UAVLocX / UAVLocY);
						state = 0;
					}
				}
				// Then let the UAVs fly along the meridians around the sphere.
				else
				{
					theta += speed * 0.01 * 1.05 / (PI * 10);
					UAVLocZ = 50 + 10 * sin(theta);
					UAVLocY = sqrt((pow(10, 2) - pow(UAVLocZ - 50, 2)) / (pow(ratioXY, 2) + 1));
					UAVLocX = UAVLocY * ratioXY;
					UAVLocX *= signX * cos(theta) / fabs(cos(theta));
					UAVLocY *= signY * cos(theta) / fabs(cos(theta));
				}
			}

			count += 1;
			clockTime += 0.01;
		}
	}
	
	return 0;
}