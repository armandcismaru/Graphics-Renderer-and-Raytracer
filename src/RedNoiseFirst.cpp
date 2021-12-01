#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <Colour.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath> 
#include <map>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <chrono>
#include <thread>
#include <math.h>

#define WIDTH 320 * 2
#define HEIGHT 240 * 2
#define SCALING_CONSTANT 200
#define THETA 0.261799388
#define SCALER 0.17

float pi = 2 * acos(0.0);
int rendering = 0;
bool orbit = false;

glm::mat3 camOrientation = glm::mat3(1.0);
//glm::vec3 camPos = glm::vec3(0.0, 0.0, 3.0);

glm::vec3 camPos = glm::vec3(0.0, 0.35, 0.7);
//glm::vec3 light = glm::vec3(0.0, 0.4, 0.00);
glm::vec3 light = glm::vec3(0.3, 0.35, 0.5);

std::array<std::array<float, HEIGHT>, WIDTH> depthBuffer;
std::vector<std::array<glm::vec3, 3>> faceNormals;

void drawStroked(DrawingWindow &window, CanvasTriangle triangle, Colour c);
void drawFilled(DrawingWindow &window, CanvasTriangle triangle, Colour c);

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues)
{
	std::vector<float> result;
	float delta = (to - from) / (numberOfValues);
	for (int i = 0; i < numberOfValues; i++)
	{
		result.push_back(from + i * delta);
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues)
{
	std::vector<glm::vec3> result;
	std::vector<float> one;
	std::vector<float> two;
	std::vector<float> three;
	one = interpolateSingleFloats(from[0], to[0], numberOfValues);
	two = interpolateSingleFloats(from[1], to[1], numberOfValues);
	three = interpolateSingleFloats(from[2], to[2], numberOfValues);

	for (int i = 0; i < numberOfValues; i++)
	{
		glm::vec3 triple = glm::vec3(one[i], two[i], three[i]);
		result.push_back(triple);
	}
	return result;
}

void drawRGB(DrawingWindow &window)
{
	glm::vec3 topLeft(255, 0, 0);	   // red
	glm::vec3 topRight(0, 0, 255);	   // blue
	glm::vec3 bottomRight(0, 255, 0);  // green
	glm::vec3 bottomLeft(255, 255, 0); // yellow

	std::vector<glm::vec3> firstCol = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT);
	std::vector<glm::vec3> lastCol = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);

	window.clearPixels();
	std::vector<glm::vec3> leftRight;

	for (size_t y = 0; y < window.height; y++)
	{
		leftRight = interpolateThreeElementValues(firstCol[y], lastCol[y], WIDTH);
		for (size_t x = 0; x < window.width; x++)
		{
			glm::vec3 pixel = leftRight[x];
			float red = pixel[0];
			float green = pixel[1];
			float blue = pixel[2];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void drawGreyscale(DrawingWindow &window)
{
	window.clearPixels();
	std::vector<float> scales;
	scales = interpolateSingleFloats(255, 0, window.height * window.width);
	int count = 0;

	for (size_t x = 0; x < window.width; x++)
	{
		for (size_t y = 0; y < window.height; y++)
		{
			float red = scales[count];
			float green = scales[count];
			float blue = scales[count];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			count++;
			window.setPixelColour(x, y, colour);
		}
	}
}

glm::mat3 xRotationMatrix(float theta)
{
	glm::mat3 mat = glm::mat3(
		glm::vec3(1.0, 0.0, 0.0),
		glm::vec3(0.0, std::cos(theta), std::sin(theta)),
		glm::vec3(0.0, -std::sin(theta), std::cos(theta)));
	return mat;
}

glm::mat3 yRotationMatrix(float theta)
{
	glm::mat3 mat = glm::mat3(
		glm::vec3(std::cos(theta), 0.0, -std::sin(theta)),
		glm::vec3(0.0, 1.0, 0.0),
		glm::vec3(std::sin(theta), 0.0, std::cos(theta)));
	return mat;
}

void handleEvent(SDL_Event event, DrawingWindow &window)
{
	if (event.type == SDL_KEYDOWN)
	{
		if (event.key.keysym.sym == SDLK_LEFT)
			camPos[0] += 0.1;
		else if (event.key.keysym.sym == SDLK_RIGHT)
			camPos[0] -= 0.1;
		else if (event.key.keysym.sym == SDLK_UP)
			camPos[1] -= 0.1;
		else if (event.key.keysym.sym == SDLK_DOWN)
			camPos[1] += 0.1;
		else if (event.key.keysym.sym == SDLK_i)
			light[2] -= 0.05;
		else if (event.key.keysym.sym == SDLK_k)
			light[2] += 0.05;
		else if (event.key.keysym.sym == SDLK_j)
			light[0] -= 0.05;
		else if (event.key.keysym.sym == SDLK_l)
			light[0] += 0.05;
		else if (event.key.keysym.sym == SDLK_m)
		{
			rendering += 1;
			rendering %= 3;
		}
		else if (event.key.keysym.sym == SDLK_u)
		{
			CanvasPoint v0 = CanvasPoint(std::rand() % WIDTH, std::rand() % HEIGHT);
			CanvasPoint v1 = CanvasPoint(std::rand() % WIDTH, std::rand() % HEIGHT);
			CanvasPoint v2 = CanvasPoint(std::rand() % WIDTH, std::rand() % HEIGHT);
			Colour col = Colour(std::rand() % 255, std::rand() % 255, std::rand() % 255);
			drawStroked(window, CanvasTriangle(v0, v1, v2), col);
		}
		else if (event.key.keysym.sym == SDLK_f)
		{
			CanvasPoint v0 = CanvasPoint(std::rand() % WIDTH, std::rand() % HEIGHT);
			CanvasPoint v1 = CanvasPoint(std::rand() % WIDTH, std::rand() % HEIGHT);
			CanvasPoint v2 = CanvasPoint(std::rand() % WIDTH, std::rand() % HEIGHT);
			Colour col = Colour(std::rand() % 255, std::rand() % 255, std::rand() % 255);
			drawFilled(window, CanvasTriangle(v0, v1, v2), col);
		}
		else if (event.key.keysym.sym == SDLK_p)
		{
			if (orbit)
				orbit = false;
			else
				orbit = true;
		}
		else if (event.key.keysym.sym == SDLK_w)
		{
			glm::mat3 matX = xRotationMatrix(THETA);
			camPos = matX * camPos;
			camOrientation = camOrientation * matX;
		}
		else if (event.key.keysym.sym == SDLK_s)
		{
			glm::mat3 matX = xRotationMatrix(-THETA);
			camPos = matX * camPos;
			camOrientation = camOrientation * matX;
		}
		else if (event.key.keysym.sym == SDLK_q)
		{
			glm::mat3 matY = yRotationMatrix(THETA);
			camPos = matY * camPos;
			camOrientation = matY * camOrientation;
		}
		else if (event.key.keysym.sym == SDLK_e)
		{
			glm::mat3 matY = yRotationMatrix(-THETA);
			camPos = matY * camPos;
			camOrientation = matY * camOrientation;
		}
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN)
	{
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	else if (event.type == SDL_MOUSEWHEEL)
	{
		if (event.wheel.y > 0)
		{
			camPos[2] -= 0.5;
		}
		else if (event.wheel.y < 0)
		{
			camPos[2] += 0.5;
		}
	}
}

void drawLine(CanvasPoint from, CanvasPoint to, Colour c, DrawingWindow &window)
{
	if ((from.x < 0 || to.x >= WIDTH) && (to.y >= HEIGHT || from.y < 0))
		return;

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	int numberOfSteps = std::max(abs(xDiff), abs(yDiff));

	std::vector<glm::vec3> values = interpolateThreeElementValues(
		glm::vec3(floor(from.x), from.y, from.depth),
		glm::vec3(ceil(to.x), to.y, to.depth),
		numberOfSteps + 1);

	for (int i = 0; i < values.size(); i++)
	{
		int x = values[i][0];
		int y = round(values[i][1]);
		float depth = values[i][2];
		float a, b;

		uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;

		if (depth < 0 && depthBuffer[x][y] < 0)
		{
			a = abs(depth);
			b = abs(depthBuffer[x][y]);
		}
		else
		{
			a = depth;
			b = depthBuffer[x][y];
		}

		if (a >= b || depthBuffer[x][y] == 0)
		{
			depthBuffer[x][y] = depth;
			if ((x >= 0 && x < WIDTH) && (y < HEIGHT && y >= 0))
				window.setPixelColour(x, y, colour);
		}
	}
}

void drawStroked(DrawingWindow &window, CanvasTriangle triangle, Colour c)
{
	drawLine(triangle.v0(), triangle.v1(), c, window);
	drawLine(triangle.v0(), triangle.v2(), c, window);
	drawLine(triangle.v1(), triangle.v2(), c, window);
}

void fillTopFlat(DrawingWindow &window, CanvasPoint v1, CanvasPoint v2, CanvasPoint v3, Colour c)
{
	float invSlope1 = (v3.x - v1.x) / (v3.y - v1.y);
	float invSlope2 = (v3.x - v2.x) / (v3.y - v1.y);
	float invSlopeD1 = (v3.depth - v1.depth) / (v3.y - v1.y);
	float invSlopeD2 = (v3.depth - v2.depth) / (v3.y - v2.y);

	float lineLeft = v3.x;
	float lineRight = v3.x;
	float depthLeft = v3.depth;
	float depthRight = v3.depth;

	for (int y = v3.y; y >= v1.y; y--)
	{
		drawLine(CanvasPoint(lineLeft, y, depthLeft), CanvasPoint(lineRight, y, depthRight), c, window);
		lineLeft -= invSlope1;
		lineRight -= invSlope2;
		depthLeft -= invSlopeD1;
		depthRight -= invSlopeD2;
	}
}

void fillBottomFlat(DrawingWindow &window, CanvasPoint v1, CanvasPoint v2, CanvasPoint v3, Colour c)
{
	float invSlope1 = (v2.x - v1.x) / (v2.y - v1.y);
	float invSlope2 = (v3.x - v1.x) / (v2.y - v1.y);
	float invSlopeD1 = (v2.depth - v1.depth) / (v2.y - v1.y);
	float invSlopeD2 = (v3.depth - v1.depth) / (v2.y - v1.y);

	float lineLeft = v1.x;
	float lineRight = v1.x;
	float depthLeft = v1.depth;
	float depthRight = v1.depth;

	for (int y = v1.y; y <= v2.y; y++)
	{
		drawLine(CanvasPoint(lineLeft, y, depthLeft), CanvasPoint(lineRight, y, depthRight), c, window);
		lineLeft += invSlope1;
		lineRight += invSlope2;
		depthLeft += invSlopeD1;
		depthRight += invSlopeD2;
	}
}

void drawFilled(DrawingWindow &window, CanvasTriangle triangle, Colour c)
{
	for (int i = 0; i < 2; i++)
		for (int j = i + 1; j < 3; j++)
			if (triangle[i].y > triangle[j].y)
				std::swap(triangle[i], triangle[j]);

	if (round(triangle.v1().y) == round(triangle.v2().y))
	{
		fillBottomFlat(window, triangle.v0(), triangle.v1(), triangle.v2(), c);
	}
	else if (round(triangle.v0().y) == round(triangle.v1().y))
	{
		fillTopFlat(window, triangle.v0(), triangle.v1(), triangle.v2(), c);
	}
	else
	{
		float v3x = triangle.v0().x + ((triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y)) * (triangle.v2().x - triangle.v0().x);

		float v3depth = triangle.v0().depth + ((triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y)) * (triangle.v2().depth - triangle.v0().depth);

		CanvasPoint v3 = CanvasPoint(v3x, triangle.v1().y, v3depth);
		fillBottomFlat(window, triangle.v0(), triangle.v1(), v3, c);
		fillTopFlat(window, triangle.v1(), v3, triangle.v2(), c);
	}
}

std::vector<ModelTriangle> readOBJ(std::string fileName, std::map<std::string, Colour> colours, DrawingWindow &window)
{
	std::ifstream file(fileName);
	std::string line;
	std::vector<ModelTriangle> triangles;
	std::vector<glm::vec3> ventrices;
	std::vector<glm::vec3> normals;

	if (file.is_open())
	{
		std::string cur;
		while (file)
		{
			std::getline(file, line);
			std::vector<std::string> c = split(line, ' ');

			if (c[0] == "usemtl")
			{
				cur = c[1];
			}
			else if (c[0] == "vn")
			{
				float x = std::stod(c[1]);
				float y = std::stod(c[2]);
				float z = std::stod(c[3]);
				normals.push_back(glm::vec3(x, y ,z));
			}
			else if (c[0] == "v")
			{
				float x = std::stod(c[1]) * SCALER;
				float y = std::stod(c[2]) * SCALER;
				float z = std::stod(c[3]) * SCALER;
				ventrices.push_back(glm::vec3(x, y, z));
			}
			else if (c[0] == "f")
			{
				int xn=0; int yn=0; int zn=0;
				Colour color = Colour();
				std::vector<std::string> ind = split(line, ' ');

				for (int i = 1; i <= 3; i++)
					ind[i].erase(std::remove(ind[i].begin(), ind[i].end(), '/'), ind[i].end());
			
				/*xn = std::stoi(ind[1]) - 1;
				yn = std::stoi(ind[2]) - 1;
				zn = std::stoi(ind[3]) - 1;*/

				xn = std::stoi(ind[1].substr(0, ind[1].size()/2)) - 1;
				yn = std::stoi(ind[2].substr(0, ind[2].size()/2)) - 1;
				zn = std::stoi(ind[3].substr(0, ind[3].size()/2)) - 1;
				std::array<glm::vec3, 3> norms = {normals[xn], normals[yn], normals[zn]};
				faceNormals.push_back(norms);

				if (!colours.empty())
					color = Colour(cur, colours[cur].red, colours[cur].green, colours[cur].blue);
				else 
					color = Colour("Red", 201, 22, 22);

				ModelTriangle t = ModelTriangle(ventrices[xn], ventrices[yn], ventrices[zn], color); 
				t.normal = glm::normalize(glm::cross(ventrices[xn] - ventrices[zn], ventrices[yn] - ventrices[zn]));
				triangles.push_back(t);
			}
		}
	}
	return triangles;
}

std::map<std::string, Colour> readMTL(std::string fileName)
{
	std::ifstream file(fileName);
	std::string line;
	std::map<std::string, Colour> colours;

	if (file.is_open())
	{
		while (file)
		{
			std::getline(file, line);
			std::vector<std::string> ind = split(line, ' ');

			if (ind[0] == "newmtl")
			{
				std::string color = ind[1];
				std::getline(file, line);
				std::vector<std::string> rgb = split(line, ' ');
				int r = stof(rgb[1]) * 255;
				int g = stof(rgb[2]) * 255;
				int b = stof(rgb[3]) * 255;
				colours[color] = Colour(color, r, g, b);
			}
		}
	}
	return colours;
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength)
{
	glm::vec3 cameraToVertex = vertexPosition - cameraPosition;
	glm::vec3 adjustedCamera = glm::normalize(cameraToVertex * camOrientation);

	float u = focalLength * (adjustedCamera[0] / adjustedCamera[2]) * SCALING_CONSTANT + WIDTH / 2;
	float v = focalLength * (adjustedCamera[1] / adjustedCamera[2]) * SCALING_CONSTANT + HEIGHT / 2;

	return CanvasPoint(WIDTH - u, v, 1 / (cameraToVertex[2]));
}

void renderTriangles(std::vector<ModelTriangle> triangles, glm::vec3 cameraPosition, float focalLength, DrawingWindow &window)
{
	for (int i = 0; i < triangles.size(); i++)
	{
		ModelTriangle triangle = triangles[i];
		std::array<CanvasPoint, 3> cPoints{};

		for (int j = 0; j < 3; j++)
		{
			cPoints[j] = getCanvasIntersectionPoint(cameraPosition, triangle.vertices[j], focalLength);
		}
		CanvasTriangle strokedTriangle = CanvasTriangle(cPoints[0], cPoints[1], cPoints[2]);
		switch (rendering)
		{
		case 1:
			drawFilled(window, strokedTriangle, triangle.colour);
			break;
		case 2:
			drawStroked(window, strokedTriangle, triangle.colour);
			break;
		}
	}
}

glm::mat3 lookAt(glm::vec3 point)
{
	glm::vec3 zaxis = glm::normalize(point - camPos);
	glm::vec3 xaxis = glm::normalize(glm::cross(zaxis, glm::vec3(0, 1, 0)));
	glm::vec3 yaxis = glm::cross(xaxis, zaxis);

	zaxis = -zaxis;

	glm::mat3 viewMatrix = glm::mat3(
		glm::vec3(xaxis.x, xaxis.y, xaxis.z),
		glm::vec3(yaxis.x, yaxis.y, yaxis.z),
		glm::vec3(zaxis.x, zaxis.y, zaxis.z));

	return viewMatrix;
}

void drawRasterised(DrawingWindow &window, std::vector<ModelTriangle> triangles)
{
	window.clearPixels();
	for (int i = 0; i < WIDTH; i++)
		for (int j = 0; j < HEIGHT; j++)
			depthBuffer[i][j] = 0.0;

	float focalLength = 1.5;
	renderTriangles(triangles, camPos, focalLength, window);
	window.renderFrame();

	if (orbit)
	{
		glm::mat3 matX = yRotationMatrix(0.0174532925);
		camPos = matX * camPos;
		camOrientation = lookAt(glm::vec3(0, 0, 0));
	}
}

std::vector<RayTriangleIntersection> getClosestIntersection(std::vector<ModelTriangle> triangles, glm::vec3 cameraPosition, glm::vec3 rayDirection)
{
	float tDist = INT32_MAX;
	std::vector<RayTriangleIntersection> closestIntersection;

	for (int i = 0; i < triangles.size(); i++)
	{
		ModelTriangle triangle = triangles[i];
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		float t = possibleSolution[0];
		float u = possibleSolution[1];
		float v = possibleSolution[2];

		if ((u >= 0.0 && u <= 1.0) && (v >= 0.0 && v <= 1.0) && (u + v <= 1.0))
		{
			if (t < tDist && t > 0.00001f)
			{
				tDist = t;
				glm::vec3 point = triangle.vertices[0] + u * e0 + v * e1;
				// glm::vec3 point = camPos + rayDirection * t;
				closestIntersection.push_back(RayTriangleIntersection(point, t, triangle, i));
			}
		}
	}
	return closestIntersection;
}

void printIntersections(std::vector<RayTriangleIntersection> inter)
{
	for (int i = 0; i < inter.size(); i++)
		std::cout << inter[i].intersectedTriangle.colour << '\n';
	std::cout << '\n';
}

bool intersectionHas(std::vector<RayTriangleIntersection> inter, std::string s)
{
	for (int i = 0; i < inter.size(); i++)
		if (inter[i].intersectedTriangle.colour.name == s)
			return true;
	return false;
}

bool isShadow(glm::vec3 point, std::vector<ModelTriangle> triangles, int index)
{
	std::vector<RayTriangleIntersection> lightPoints = getClosestIntersection(triangles, point, glm::normalize(light - point));
	// If light too far away it will not intersect white but still cast shadow
	if (!lightPoints.empty() && intersectionHas(lightPoints, "White"))
	{
		int shadow_index = lightPoints.back().triangleIndex;
		if (index == shadow_index && lightPoints[lightPoints.size() - 2].intersectedTriangle.colour.name != "White")
			return true;
		if (lightPoints.back().intersectedTriangle.colour.name != "White" && index != shadow_index)
			return true;
	}
	return false;
}

float getBrightness(glm::vec3 point)
{
	float length = glm::length(point - light);
	return 1 / (4 * pi * length * length);
}

glm::vec3 barycentric(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 point)
{
	glm::vec3 vec1 = v2 - v1;
	glm::vec3 vec2 = v3 - v1;
	glm::vec3 vec3 = point - v1;

	float d00 = glm::dot(vec1, vec1);
	float d01 = glm::dot(vec1, vec2);
	float d11 = glm::dot(vec2, vec2);
	float d20 = glm::dot(vec3, vec1);
	float d21 = glm::dot(vec3, vec2);
	float denom = d00 * d11 - d01 * d01;
	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	float u = 1.0f - v - w;

	return glm::vec3(u, v, w);
}

float normalBrightness (glm::vec3 point, glm::vec3 normal, bool isShadow)
{
	float brightness = getBrightness(point);
	float angle = abs(glm::dot(normal, glm::normalize(light-point)));
	
	float coef = 2 * glm::dot(glm::normalize(light-point), normal);
	glm::vec3 reflection = glm::normalize(light-point) - coef * normal;
	float spread = std::pow(glm::dot(glm::normalize(camPos-point), reflection), 256) ;

	float multiplier = brightness * angle + brightness * abs(spread) + 0.2f;
	//float multiplier = std::max(brightness * angle + brightness * abs(spread), 0.2f);
	if (isShadow) multiplier = brightness * angle * 0.5;

	return multiplier;
}

float gouraudBrightness(glm::vec3 point, ModelTriangle &triangle, int index)
{
	glm::vec3 normal1 = faceNormals[index][0];
	glm::vec3 normal2 = faceNormals[index][1];
	glm::vec3 normal3 = faceNormals[index][2];

	glm::vec3 bary = barycentric(triangle.vertices[0], triangle.vertices[1], triangle.vertices[2], point);
	
	float b1 = getBrightness(triangle.vertices[0]);
	float b2 = getBrightness(triangle.vertices[1]);
	float b3 = getBrightness(triangle.vertices[2]);

	float angle1 = abs(glm::dot(normal1, glm::normalize(light-triangle.vertices[0]))); 
	float angle2 = abs(glm::dot(normal2, glm::normalize(light-triangle.vertices[1]))); 
	float angle3 = abs(glm::dot(normal3, glm::normalize(light-triangle.vertices[2]))); 

	float coef1 = 2 * glm::dot(glm::normalize(light-triangle.vertices[0]), normal1);
	float coef2 = 2 * glm::dot(glm::normalize(light-triangle.vertices[1]), normal2);
	float coef3 = 2 * glm::dot(glm::normalize(light-triangle.vertices[2]), normal3);

	glm::vec3 reflection1 = glm::normalize(light-triangle.vertices[0]) - coef1 * normal1;
	glm::vec3 reflection2 = glm::normalize(light-triangle.vertices[1]) - coef2 * normal2;
	glm::vec3 reflection3 = glm::normalize(light-triangle.vertices[2]) - coef3 * normal3;

	float spread1 = std::pow(glm::dot(glm::normalize(camPos-triangle.vertices[0]), reflection1), 256);
	float spread2 = std::pow(glm::dot(glm::normalize(camPos-triangle.vertices[1]), reflection2), 256);
	float spread3 = std::pow(glm::dot(glm::normalize(camPos-triangle.vertices[2]), reflection3), 256);

	float m1 = b1 * angle1 + b1 * abs(spread1) + 0.3f;
	float m2 = b2 * angle2 + b2 * abs(spread2) + 0.3f;
	float m3 = b3 * angle3 + b3 * abs(spread3) + 0.3f;

	return bary[0] * m1 + bary[1] * m2 + bary[2] * m3;
}

void drawRayTraced(DrawingWindow &window, std::vector<ModelTriangle> triangles)
{
	window.clearPixels();
	float focalLength = 1.5;

	for (int i = 0; i < WIDTH; i++)
		for (int j = 0; j < HEIGHT; j++)
		{
			float x = ((i - WIDTH / 2) / (SCALING_CONSTANT * 0.8));
			float y = (((HEIGHT - j) - HEIGHT / 2) / (SCALING_CONSTANT * 0.8));

			std::vector<RayTriangleIntersection> closestPoint = getClosestIntersection(
				triangles, 
				camPos, 
				glm::normalize(glm::vec3(x, y, -focalLength) - camPos));
	
			if (!closestPoint.empty())
			{
				glm::vec3 point = closestPoint.back().intersectionPoint;
				ModelTriangle triangle = closestPoint.back().intersectedTriangle;
				int t = closestPoint.back().triangleIndex;

				float m = std::min(gouraudBrightness(point, triangle, t), 1.0f);
				//float m = std::min(normalBrightness(point, triangles[t].normal, isShadow(point, triangles, t)), 1.0f);

				Colour c = triangles[t].colour;
				uint32_t colour = (255 << 24) +
					(int(c.red * m) % 256 << 16) +
					(int(c.green * m) % 256 << 8) +
					int(c.blue * m) % 256;

				window.setPixelColour(i, j, colour);
			}
		}
	window.renderFrame();

	if (orbit)
	{
		glm::mat3 matX = yRotationMatrix(0.0174532925);
		camPos = matX * camPos;
		camOrientation = lookAt(glm::vec3(0, 0, 0));
	}
}

int main(int argc, char *argv[])
{
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::srand(std::time(nullptr));

	//std::map<std::string, Colour> colours = readMTL("cornell-box.mtl");
	std::map<std::string, Colour> colours;
	std::vector<ModelTriangle> triangles = readOBJ("sphere.obj", colours, window);

	while (true)
	{
		if (window.pollForInputEvents(event))
			handleEvent(event, window);

		switch (rendering)
		{
		case 0:
			drawRayTraced(window, triangles);
			break;
		case 1:
			drawRasterised(window, triangles);
			break;
		case 2:
			drawRasterised(window, triangles);
			break;
		}
	}
}
