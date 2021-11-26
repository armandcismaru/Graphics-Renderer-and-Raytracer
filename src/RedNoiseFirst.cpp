#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <Colour.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <chrono>
#include <thread>


#define WIDTH 320*3
#define HEIGHT 240*3
#define SCALING_CONSTANT 500
#define THETA 0.261799388

glm::vec3 camPos = glm::vec3(0.0, 0.0, 2.5);
glm::mat3 camOrientation = glm::mat3(1.0);
std::array<std::array<float, HEIGHT>, WIDTH> depthBuffer;

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues);
void drawStroked (DrawingWindow &window, CanvasTriangle triangle, Colour c);
void drawFilled (DrawingWindow &window, CanvasTriangle triangle, Colour c);
CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength);
void fillTriangleDepth(DrawingWindow &window, CanvasTriangle triangle, std::array<std::array<float, HEIGHT>, WIDTH> &depthBuffer);

void drawGreyscale(DrawingWindow &window) {
	window.clearPixels();
	std::vector<float> scales;
	scales = interpolateSingleFloats(255, 0, window.height * window.width);
	int count = 0;

	for (size_t x = 0; x < window.width; x++) {
		for (size_t y = 0; y < window.height; y++) {
			float red = scales[count];
			float green = scales[count];
			float blue = scales[count];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			count++;
			window.setPixelColour(x, y, colour);
		}
	}
}

glm::mat3 xRotationMatrix (float theta){
	glm::mat3 mat = glm::mat3(
		glm::vec3(1.0, 0.0, 0.0),
		glm::vec3(0.0, std::cos(theta), std::sin(theta)),
		glm::vec3(0.0, -std::sin(theta), std::cos(theta)) 
	);
	return mat;
}

glm::mat3 yRotationMatrix (float theta){
	glm::mat3 mat = glm::mat3(
		glm::vec3(std::cos(theta), 0.0, -std::sin(theta)),
		glm::vec3(0.0, 1.0, 0.0),
		glm::vec3(std::sin(theta), 0.0, std::cos(theta))
	);
	return mat;
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) camPos[0] += 0.1;
		else if (event.key.keysym.sym == SDLK_RIGHT) camPos[0] -= 0.1;
		else if (event.key.keysym.sym == SDLK_UP) camPos[1] -= 0.1;
		else if (event.key.keysym.sym == SDLK_DOWN) camPos[1] += 0.1;
		else if (event.key.keysym.sym == SDLK_u){
			CanvasPoint v0 = CanvasPoint(std::rand()%WIDTH, std::rand()%HEIGHT);
			CanvasPoint v1 = CanvasPoint(std::rand()%WIDTH, std::rand()%HEIGHT);
			CanvasPoint v2 = CanvasPoint(std::rand()%WIDTH, std::rand()%HEIGHT);
			Colour col = Colour(std::rand()%255, std::rand()%255, std::rand()%255);
			drawStroked(window, CanvasTriangle(v0, v1, v2), col);
		}
		else if (event.key.keysym.sym == SDLK_f){
			CanvasPoint v0 = CanvasPoint(std::rand()%WIDTH, std::rand()%HEIGHT);
			CanvasPoint v1 = CanvasPoint(std::rand()%WIDTH, std::rand()%HEIGHT);
			CanvasPoint v2 = CanvasPoint(std::rand()%WIDTH, std::rand()%HEIGHT);
			Colour col = Colour(std::rand()%255, std::rand()%255, std::rand()%255);
			drawFilled(window, CanvasTriangle(v0, v1, v2), col);
		}
		else if (event.key.keysym.sym == SDLK_w){
			glm::mat3 matX = xRotationMatrix(THETA);
			camPos = matX*camPos;
			camOrientation = camOrientation * matX;
		}
		else if (event.key.keysym.sym == SDLK_s){
			glm::mat3 matX = xRotationMatrix(-THETA);
			camPos = matX*camPos;
			camOrientation = camOrientation * matX;
		}
		else if (event.key.keysym.sym == SDLK_q){
			glm::mat3 matY = yRotationMatrix(THETA);
			camPos = matY*camPos;
			camOrientation = camOrientation * matY;
		}
		else if (event.key.keysym.sym == SDLK_e){
			glm::mat3 matY = yRotationMatrix(-THETA);
			camPos = matY*camPos;
			camOrientation = camOrientation * matY;
		}
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	} 
	else if (event.type == SDL_MOUSEWHEEL){
		if (event.wheel.y > 0){
			camPos[2] -= 0.5;
		}
		else if (event.wheel.y < 0){
			camPos[2] += 0.5;
		}
	}
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues){
	std::vector<float> result;
	float delta = (to-from)/(numberOfValues);
    for(int i=0; i<numberOfValues; i++) {
        result.push_back(from + i*delta);
    }
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues){
	std::vector<glm::vec3> result;
	std::vector<float> one;
	std::vector<float> two;
	std::vector<float> three;
	one = interpolateSingleFloats(from[0], to[0], numberOfValues);
	two = interpolateSingleFloats(from[1], to[1], numberOfValues);
	three = interpolateSingleFloats(from[2], to[2], numberOfValues);
	
	for(int i=0; i<numberOfValues; i++) {
		glm::vec3 triple = glm::vec3(one[i], two[i], three[i]);
        result.push_back(triple);
    }
	return result;
}

void drawRGB(DrawingWindow &window) {
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow

	std::vector<glm::vec3> firstCol = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT);
	std::vector<glm::vec3> lastCol = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);

	window.clearPixels();
	std::vector<glm::vec3> leftRight; 

	for (size_t y = 0; y < window.height; y++) {
		leftRight = interpolateThreeElementValues(firstCol[y], lastCol[y], WIDTH);
		for (size_t x = 0; x < window.width; x++) {
			glm::vec3 pixel = leftRight[x];
			float red = pixel[0];
			float green = pixel[1];
			float blue = pixel[2];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void drawLine(CanvasPoint from, CanvasPoint to, Colour c, DrawingWindow &window){
	if ((from.x < 0 || to.x >= WIDTH) && (to.y >= HEIGHT || from.y < 0))
		return;

	float xDiff = to.x - from.x;
	//float yDiff = to.y - from.y;
	//int numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	int numberOfSteps = abs(xDiff);
	std::vector<glm::vec3> values = interpolateThreeElementValues(
		glm::vec3(floor(from.x), from.y, from.depth), 
		glm::vec3(ceil(to.x), to.y, to.depth), 
		numberOfSteps+1
	);

	for(int i=0; i<numberOfSteps+1; i++){
		int x = values[i][0];
		int y = round(values[i][1]);
		float depth = values[i][2];
		float a, b;

		uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;

		if (depth < 0 && depthBuffer[x][y] < 0) {
			a = abs(depth);
			b = abs(depthBuffer[x][y]);
		}
		else {
			a = depth;
			b = depthBuffer[x][y];
		}

		if (a >= b || depthBuffer[x][y] == 0){
			depthBuffer[x][y] = depth;
			if ((x >= 0 && x < WIDTH) && (y < HEIGHT && y >= 0))
				window.setPixelColour(x, y, colour);
		}
	}
}

void drawStroked(DrawingWindow &window, CanvasTriangle triangle, Colour c){
	drawLine(triangle.v0(), triangle.v1(), c, window);
	drawLine(triangle.v0(), triangle.v2(), c, window);
	drawLine(triangle.v1(), triangle.v2(), c, window);
}

void fillTopFlat(DrawingWindow &window, CanvasPoint v1, CanvasPoint v2, CanvasPoint v3, Colour c){
	float invSlope1 = (v3.x - v1.x)/(v3.y - v1.y);
	float invSlope2 = (v3.x - v2.x)/(v3.y - v1.y);
	float invSlopeD1 = (v3.depth - v1.depth)/(v3.y - v1.y);
	float invSlopeD2 = (v3.depth - v2.depth)/(v3.y - v2.y);

	float lineLeft = v3.x;
	float lineRight = v3.x;
	float depthLeft = v3.depth;
	float depthRight = v3.depth;

	for (int y = v3.y; y >= v1.y; y--){
		drawLine(CanvasPoint(lineLeft, y, depthLeft), CanvasPoint(lineRight, y, depthRight), c, window);
		lineLeft -= invSlope1;
		lineRight -= invSlope2;
		depthLeft -= invSlopeD1;
		depthRight -= invSlopeD2;
	}
}

void fillBottomFlat(DrawingWindow &window, CanvasPoint v1, CanvasPoint v2, CanvasPoint v3, Colour c){
	float invSlope1 = (v2.x - v1.x)/(v2.y - v1.y);
	float invSlope2 = (v3.x - v1.x)/(v2.y - v1.y);
	float invSlopeD1 = (v2.depth - v1.depth)/(v2.y - v1.y);
	float invSlopeD2 = (v3.depth - v1.depth)/(v2.y - v1.y);

	float lineLeft = v1.x;
	float lineRight = v1.x;
	float depthLeft = v1.depth;
	float depthRight = v1.depth;

	for (int y = v1.y; y < v2.y; y++){
		drawLine(CanvasPoint(lineLeft, y, depthLeft), CanvasPoint(lineRight, y, depthRight), c, window);
		lineLeft += invSlope1;
		lineRight += invSlope2;
		depthLeft += invSlopeD1;
		depthRight += invSlopeD2;
	}
}

void drawFilled(DrawingWindow &window, CanvasTriangle triangle, Colour c){
	for (int i=0; i<2; i++)
		for(int j=i+1; j<3; j++)
			if(triangle[i].y > triangle[j].y)
				std::swap(triangle[i], triangle[j]);
	
	if (round(triangle.v1().y) == round(triangle.v2().y)){
		fillBottomFlat(window, triangle.v0(), triangle.v1(), triangle.v2(), c);
	}
	else if (round(triangle.v0().y) == round(triangle.v1().y)){
		fillTopFlat(window, triangle.v0(), triangle.v1(), triangle.v2(), c);
	}
	else
	{
		float v3x = triangle.v0().x 
			+ ((triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y)) 
				* (triangle.v2().x - triangle.v0().x);
		
		float v3depth = triangle.v0().depth 
			+ ((triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y)) 
				* (triangle.v2().depth - triangle.v0().depth);

		CanvasPoint v3 = CanvasPoint(v3x ,triangle.v1().y, v3depth);
		fillBottomFlat(window, triangle.v0(), triangle.v1(), v3, c);
		fillTopFlat(window, triangle.v1(), v3, triangle.v2(), c);
	}
}

std::vector<ModelTriangle> readOBJ (std::string fileName, std::map<std::string, Colour> colours, DrawingWindow &window){
	std::ifstream file (fileName);
	std::string line;
	std::vector<ModelTriangle> triangles;
	std::vector<glm::vec3> ventrices;

	if (file.is_open()) {
		std::string cur;
		while (file) {
			std::getline (file, line);
			std::vector<std::string> c = split(line, ' ');

			if (c[0] == "usemtl"){
				cur = c[1];
			}
			else if (c[0] == "v"){
				float x = std::stod(c[1]) * 0.17;
				float y = std::stod(c[2]) * 0.17;
				float z = std::stod(c[3]) * 0.17;
				ventrices.push_back(glm::vec3(x, y ,z));
			}
			else if (c[0] == "f"){
				std::vector<std::string> ind = split(line, ' ');
				for (int i=1; i<=3; i++)
					ind[i].erase(std::remove(ind[i].begin(), ind[i].end(), '/'), ind[i].end());

				int xn = std::stoi(ind[1])-1;
				int yn = std::stoi(ind[2])-1;
				int zn = std::stoi(ind[3])-1;
				Colour color = Colour(cur, colours[cur].red, colours[cur].green, colours[cur].blue);
				triangles.push_back(ModelTriangle(ventrices[xn], ventrices[yn], ventrices[zn], color));
			}
		}
	}
	return triangles;
}

std::map<std::string, Colour> readMTL(std::string fileName){
	std::ifstream file(fileName);
	std::string line;
	std::map<std::string, Colour> colours;

	if (file.is_open()){
		while (file){
			std::getline(file, line);
			std::vector<std::string> ind = split(line, ' ');

			if (ind[0] == "newmtl"){
				std::string color = ind[1];
				std::getline(file, line);
				std::vector<std::string> rgb = split(line, ' ');
				int r = stof(rgb[1])*255;
				int g = stof(rgb[2])*255;
				int b = stof(rgb[3])*255;
				colours[color] = Colour(color, r, g, b);
			}	
		}
	}
	return colours;
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength){
	glm::vec3 cameraToVertex = vertexPosition - cameraPosition;
	glm::vec3 adjustedCamera = glm::normalize(cameraToVertex * camOrientation);
	
	float u = focalLength * (adjustedCamera[0]/adjustedCamera[2]) * SCALING_CONSTANT + WIDTH/2;  
	float v = focalLength * (adjustedCamera[1]/adjustedCamera[2]) * SCALING_CONSTANT + HEIGHT/2;  

	return CanvasPoint(WIDTH-u, v, 1/(cameraToVertex[2]));
}

void renderTriangles(std::vector<ModelTriangle> triangles, glm::vec3 cameraPosition, float focalLength, DrawingWindow &window){
	for (int i=0; i<triangles.size(); i++){
		ModelTriangle triangle = triangles[i];
		std::array<CanvasPoint, 3> cPoints{}; 

		for (int j=0; j<3; j++){
			cPoints[j] = getCanvasIntersectionPoint(cameraPosition, triangle.vertices[j], focalLength);
		}
		CanvasTriangle strokedTriangle = CanvasTriangle(cPoints[0], cPoints[1], cPoints[2]);
		drawFilled(window, strokedTriangle, triangle.colour);
		//drawStroked(window, strokedTriangle, Colour(0,0,0));
	}
}

void lookAt(glm::vec3 point){

	return;
}

void draw(DrawingWindow &window, std::vector<ModelTriangle> triangles) {
	window.clearPixels();
	for (int i=0; i <WIDTH; i++)
		for (int j=0; j <HEIGHT; j++)
			depthBuffer[i][j] = 0.0;

	float focalLength = camPos[2]/2;
	renderTriangles(triangles, camPos, focalLength, window);
	window.renderFrame();
	//glm::mat3 matX = yRotationMatrix(0.00174532925);
	//camPos = matX*camPos;
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::srand(std::time(nullptr));

	std::map<std::string, Colour> colours = readMTL("cornell-box.mtl");
	std::vector<ModelTriangle> triangles = readOBJ("cornell-box.obj", colours, window);

	while (true) {
		if (window.pollForInputEvents(event))
			handleEvent(event, window);
		draw(window, triangles);
	}
}
