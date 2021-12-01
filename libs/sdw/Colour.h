#pragma once

#include <iostream>

struct Colour {
	std::string name;
	int red{};
	int green{};
	int blue{};
	Colour();
	Colour(int r, int g, int b);
	Colour(std::string n, int r, int g, int b);
};

std::ostream &operator<<(std::ostream &os, const Colour &colour);

//std::array<glm::vec3, 3> norms = {normals[xn], normals[yn], normals[zn]};
// faceNormals.push_back(norms);
// std::vector<std::array<glm::vec3, 3>> faceNormals;
