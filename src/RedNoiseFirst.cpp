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


#define WIDTH 320*2
#define HEIGHT 240*2
#define SCALING_CONSTANT 600
std::array<std::array<float, HEIGHT>, WIDTH> depthBuffer;

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues);
void drawStroked (DrawingWindow &window, CanvasTriangle triangle, Colour c);
void drawFilled (DrawingWindow &window, CanvasTriangle triangle, Colour c);
CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength);
void fillTriangleDepth(DrawingWindow &window, CanvasTriangle triangle, std::array<std::array<float, HEIGHT>, WIDTH> &depthBuffer);


void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 64.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

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

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
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

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues){
	std::vector<float> result;
	float delta = (to-from)/(numberOfValues*1.0);
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
	int count = 0;
	for(int i=0; i<numberOfValues; i++) {
		glm::vec3 triple = glm::vec3(one[count], two[count], three[count]);
        result.push_back(triple);
		count++;
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
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));

	std::vector<glm::vec3> values = interpolateThreeElementValues(
		glm::vec3(from.x, from.y, from.depth), 
		glm::vec3(to.x, to.y, to.depth), numberOfSteps+1
	);
	
	for(int i=0; i<numberOfSteps; i++){
		int x = round(values[i][0]);
		int y = round(values[i][1]);
		float depth = values[i][2];

		uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;
	
		if (depth > depthBuffer[x][y] || depthBuffer[x][y] == 0){
			depthBuffer[x][y] = depth;
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
	float invSlope2 = (v3.x - v2.x)/(v3.y - v2.y);

	float lineLeft = v3.x;
	float lineRight = v3.x;
	//std::cout << v1.y << " " << v3.y << std::endl;

	std::vector<glm::vec3> values = interpolateThreeElementValues(
		glm::vec3(v1.depth, v2.depth, v1.y), 
		glm::vec3(v3.depth, v3.depth, v3.y), 
		v3.y - v1.y + 1
	); 

	for (int y = values.size(); y >= 0; y--){
		glm::vec3 val = values[y];
		//std::cout << lDepths[y] << " " << rDepths[y] << std::endl;
		drawLine(CanvasPoint(lineLeft, val[2], val[0]), CanvasPoint(lineRight, val[2], val[1]), c, window);
		lineLeft -= invSlope1;
		lineRight -= invSlope2;
	}
}

void fillBottomFlat(DrawingWindow &window, CanvasPoint v1, CanvasPoint v2, CanvasPoint v3, Colour c){
	float invSlope1 = (v2.x - v1.x)/(v2.y - v1.y);
	float invSlope2 = (v3.x - v1.x)/(v3.y - v1.y);

	float lineLeft = v1.x;
	float lineRight = v1.x;

	std::vector<glm::vec3> values = interpolateThreeElementValues(
		glm::vec3(v1.depth, v1.depth, v1.y), 
		glm::vec3(v2.depth, v3.depth, v2.y), v2.y - v1.y + 1
	); 

	for (int y = 0; y < values.size(); y++){
		glm::vec3 val = values[y];
		drawLine(CanvasPoint(lineLeft, val[2], val[0]), CanvasPoint(lineRight, val[2], val[1]), c, window);
		lineLeft += invSlope1;
		lineRight += invSlope2;
	}
}

void drawFilled(DrawingWindow &window, CanvasTriangle triangle, Colour c){
	for (int i=0; i<2; i++)
		for(int j=i+1; j<3; j++)
			if(triangle[i].y > triangle[j].y)
				std::swap(triangle[i], triangle[j]);

	if (triangle.v1().y == triangle.v2().y){
		fillBottomFlat(window, triangle.v0(), triangle.v1(), triangle.v2(), c);
	}
	else if (triangle.v0().y == triangle.v1().y){
		fillTopFlat(window, triangle.v0(), triangle.v1(), triangle.v2(), c);
	}
	else
	{
		float v3x = triangle.v0().x 
			+ ((triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y)) 
				* (triangle.v2().x - triangle.v0().x);

		std::vector<float> leftDepths = interpolateSingleFloats(
			triangle.v0().depth, 
			triangle.v2().depth, 
			triangle.v2().y-triangle.v0().y+1
		);
		float v3depth = leftDepths[round(triangle.v1().y-triangle.v0().y+1)];

		CanvasPoint v3 = CanvasPoint(v3x ,triangle.v1().y, v3depth);
		fillBottomFlat(window, triangle.v0(), triangle.v1(), v3, c);
		fillTopFlat(window, triangle.v1(), v3, triangle.v2(), c);
	}
	//drawStroked(window, triangle, Colour(255, 255, 255));
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
				float z = std::stod(c[3]) * 0.40;
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
	float cameraX = vertexPosition[0] - cameraPosition[0];
	float cameraY = vertexPosition[1] - cameraPosition[1];
	float cameraZ = vertexPosition[2] - cameraPosition[2];

	float u = focalLength * (cameraX/cameraZ) * SCALING_CONSTANT + WIDTH/2;  
	float v = focalLength * (cameraY/cameraZ) * SCALING_CONSTANT + HEIGHT/2;  

	return CanvasPoint(WIDTH-u, v, 1/vertexPosition[2]);
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
		//drawStroked(window, strokedTriangle, triangle.colour);
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::srand(std::time(nullptr));

	std::map<std::string, Colour> colours = readMTL("cornell-box.mtl");
	std::vector<ModelTriangle> triangles = readOBJ("cornell-box.obj", colours, window);

	for (int i=0; i <WIDTH; i++)
		for (int j=0; j <HEIGHT; j++)
			depthBuffer[i][j] = 0.0;

	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.15);
	float focalLength = 2.3;

	renderTriangles(triangles, cameraPosition, focalLength, window);

	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		window.renderFrame();
		//std::this_thread::sleep_for(std::chrono::milliseconds(200));
	}
}
