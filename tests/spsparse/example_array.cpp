#include <array>
#include <iostream>

int main(int argc, char **argv)
{
	std::array<int, 3> xxx = {{5,6,2}};
	for (int i=0; i<xxx.size(); ++i) {
		printf("x[%d] = %d\n", i, xxx[i]);
	}
}
