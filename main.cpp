#include <string>
#include <iostream>
#include <tetgen.h>
#include <Eigen/Dense>
#include <unordered_map>
#include <igl/readMESH.h>
#include <igl/writeMSH.h>

void outputToMSH(const std::string& meshFile, const std::string& mshFile)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi T, F;
	igl::readMESH(meshFile, V, T, F);
}

void outputTetCoord(const tetgenio& out, int tetIdx)
{
	std::array<int, 4> cornerIndex = {
		out.tetrahedronlist[tetIdx * 4] ,
		out.tetrahedronlist[tetIdx * 4 + 1],
		out.tetrahedronlist[tetIdx * 4 + 2],
		out.tetrahedronlist[tetIdx * 4 + 3]
	};

	for (int j = 0; j < 4; ++j)
	{
		std::cout << cornerIndex[j] << " {" <<
			out.pointlist[cornerIndex[j] * 3] << ", " <<
			out.pointlist[cornerIndex[j] * 3 + 1] << ", " <<
			out.pointlist[cornerIndex[j] * 3 + 2] << "}\n";
	}
}

void outputTetFace(const tetgenio& out, int tetIdx)
{
	for (int j = 0; j < 4; ++j)
	{
		std::cout << " {" <<
			out.trifacelist[tetIdx * 12 + j * 3] << ", " <<
			out.trifacelist[tetIdx * 12 + j * 3 + 1] << ", " <<
			out.trifacelist[tetIdx * 12 + j * 3 + 2] << "}\n";
	}
}

int getNeighborTetType(const tetgenio& out, int tetIdx, int tetNeighborType)
{
	std::array<int, 4> cornerIndex = {
		out.tetrahedronlist[tetIdx * 4] ,
		out.tetrahedronlist[tetIdx * 4 + 1],
		out.tetrahedronlist[tetIdx * 4 + 2],
		out.tetrahedronlist[tetIdx * 4 + 3]
	};

	int neighborTetIdx = out.neighborlist[tetIdx * 4 + tetNeighborType];
	if (neighborTetIdx == -1)
	{
		std::cout << "neighbor #" << tetNeighborType << " is boundary neighbor.\n";
		return -1;
	}
	else
	{
		std::unordered_map<int, int> cornerToIndex;
		for (int i = 0; i < 4; ++i)
			if (i != tetNeighborType) cornerToIndex[cornerIndex[i]] = i;

		std::cout << "neighbor #" << tetNeighborType << " tet coordinate: \n";
		outputTetCoord(out, neighborTetIdx);
		std::array<int, 4> neighborTetCornerIdx = {
			out.tetrahedronlist[neighborTetIdx * 4] ,
			out.tetrahedronlist[neighborTetIdx * 4 + 1],
			out.tetrahedronlist[neighborTetIdx * 4 + 2],
			out.tetrahedronlist[neighborTetIdx * 4 + 3]
		};
		const auto predIter = std::find_if(neighborTetCornerIdx.begin(), neighborTetCornerIdx.end(), [&](int idx) {
			return std::find(cornerIndex.begin(), cornerIndex.end(), idx) == cornerIndex.end();
			});
		int typeInNeighborTet = predIter - neighborTetCornerIdx.begin();

		for (int i = 0; i < 4; ++i)
			if (i != typeInNeighborTet)
				std::cout << "neighbor tet's corner idx: " << neighborTetCornerIdx[i] <<
				" in tet's corner idx:" << cornerToIndex[neighborTetCornerIdx[i]] << "\n";

		return typeInNeighborTet;
	}
}

bool outputOneTet(const std::array<Eigen::Vector3d, 4>& tet, const int tetIdx, std::ofstream& out, int& id) {
	out << "# tetIdx = " << tetIdx << "\n";

	out << "v " << tet[0].transpose() << "\n";
	out << "v " << tet[1].transpose() << "\n";
	out << "v " << tet[2].transpose() << "\n";
	out << "v " << tet[3].transpose() << "\n";
	out << "f " << id << " " << id + 1 << " " << id + 2 << "\n";
	out << "f " << id << " " << id + 1 << " " << id + 3 << "\n";
	out << "f " << id << " " << id + 2 << " " << id + 3 << "\n";
	out << "f " << id + 1 << " " << id + 2 << " " << id + 3 << "\n";

	id += 4;

	return true;
}

bool outputTets(const tetgenio& tetgenOut, const std::string& filename) {
	std::ofstream out(filename);

	int id = 1;
	for (int i = 0; i < tetgenOut.numberoftetrahedra; ++i) {
		std::array<int, 4> cornerIndex = {
				tetgenOut.tetrahedronlist[i * 4], tetgenOut.tetrahedronlist[i * 4 + 1],
				tetgenOut.tetrahedronlist[i * 4 + 2], tetgenOut.tetrahedronlist[i * 4 + 3]
		};

		std::array<Eigen::Vector3d, 4> tet;
		for (int j = 0; j < 4; ++j)
			tet[j] = Eigen::Vector3d(tetgenOut.pointlist[cornerIndex[j] * 3],
				tetgenOut.pointlist[cornerIndex[j] * 3 + 1],
				tetgenOut.pointlist[cornerIndex[j] * 3 + 2]);

		outputOneTet(tet, i, out, id);
	}
	out.close();

	return true;
}

int main()
{
	//char paras[] = "pgnq1.414a1e-6"; // p就是最普通的四面体化，n代表输出neighbor，q代表外接球半径与最长边比，a代表体积
	char paras[] = "pgn"; // 
	tetgenbehavior arguments;

	arguments.parse_commandline(paras);

	tetgenio in, out;
	//in.load_off(const_cast<char*>("bunny_box.off"));

	tetgenmesh m;
	m.in = &in;
	m.b = &arguments;
	m.b->no_sort = 1;
	m.b->verbose = 1;

	//in.load_off(const_cast<char*>("cube.off"));
	m.in->load_ply(const_cast<char*>("cube.ply"));
	clock_t ck;
	m.incrementaldelaunay(ck);
	m.outelements(&out);
	m.outnodes(&out);

	std::cout << out.numberoftetrahedra << std::endl;
	//std::cout << out.pointlist << std::endl;
	/*for (int i = 0; i < out.numberoftetrahedra; ++i)
		outputTetCoord(out, i);*/
	const std::string test_file = R"(E:\VSProjects\TetGenExample\test.obj)";
	outputTets(out, test_file);


	//tetrahedralize(&arguments, &in, &out);
	//tetrahedralize(&arguments, &in, nullptr);

	/*for (int i = 0; i < out.numberoftetrahedra; ++i)
	{
		std::cout << "tet coordinate: \n";
		outputTetCoord(out, i);
		std::cout << "\n";

		std::cout << "tri face index: \n";
		outputTetFace(out, i);
		std::cout << "\n";

		for (int j = 0; j < 4; ++j)
		{
			int neighborType = getNeighborTetType(out, i, j);
			std::cout << "neighbor type = " << neighborType << "\n\n";
		}
		std::cout << "=========\n";
		system("pause");
	}*/
	//std::cout << out.numberoftetrahedra << std::endl;
	//std::cout << out.numberofpoints << std::endl;

	return 0;
}