#include <array>
#include <yaml-cpp/yaml.h>

#include "problem/SWE/swe_definitions.hpp"

#include "bathymetry.hpp"

enum Pattern : uchar { simple = 0, zigzag = 1, checker = 2 };

struct Node {
    uint ID;
    std::array<double, 3> coord;
};

struct Element {
    uint ID;
    uint type;
    std::vector<uint> nodes;
};

struct Boundary {
    uchar type;
    std::vector<uint> nodes;

    // Only needed for tidal and flow boundaries
    double frequency;
    double forcing_factor;
    double equilibrium_argument;
};

std::ostream& operator<<(std::ostream& os, const Boundary& bd) {
    switch (bd.type) {
        case SWE::BoundaryTypes::land:
            os << "    type: land\n";
            break;
        case SWE::BoundaryTypes::tide:
            os << "    type: tide\n";
            break;
        case SWE::BoundaryTypes::flow:
            os << "    type: flow\n";
            break;
        case SWE::BoundaryTypes::function:
            os << "    type: function\n";
            break;
        case SWE::BoundaryTypes::outflow:
            os << "    type: outflow\n";
            break;
    }

    if (bd.type == SWE::BoundaryTypes::tide || bd.type == SWE::BoundaryTypes::flow) {
        os << "    frequency: " << bd.frequency << '\n'
           << "    forcing factor: " << bd.forcing_factor << '\n'
           << "    equilibrium argument: " << bd.equilibrium_argument << '\n';
    }
    return os;
}

struct MeshGeneratorInput {
    std::string mesh_name;

    /*
        p4-------p3
         |        |
         |        |
        p1-------p2
    */

    Point<2> p1;
    Point<2> p2;
    Point<2> p3;
    Point<2> p4;

    uint num_x_subdivisions;
    uint num_y_subdivisions;

    Pattern pattern;

    std::vector<Boundary> boundary;

    MeshGeneratorInput(const std::string& input_string);

    void Summarize();
};

void simple_pattern_tri(uint, uint, std::vector<Element>&);
void zigzag_pattern_tri(uint, uint, std::vector<Element>&);
void checker_pattern_tri(uint, uint, std::vector<Element>&);

void WriteBCISFile(const std::string& mesh_name, const std::vector<Boundary>& boundaries);

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage:\n"
                  << "    quad_mesh_generator <input file>\n";
        exit(1);
    }

    const MeshGeneratorInput input(std::string{argv[1]});

    uint m = input.num_x_subdivisions;
    uint n = input.num_y_subdivisions;

    double dxi  = 2.0 / m;
    double deta = 2.0 / n;

    double N1, N2, N3, N4;
    double xi, eta;

    std::vector<Node> nodes((m + 1) * (n + 1));

    for (uint i = 0; i <= m; ++i) {
        for (uint j = 0; j <= n; ++j) {
            xi  = -1.0 + dxi * i;
            eta = -1.0 + deta * j;

            N1 = (1.0 - xi) * (1.0 - eta) / 4.0;
            N2 = (1.0 + xi) * (1.0 - eta) / 4.0;
            N3 = (1.0 + xi) * (1.0 + eta) / 4.0;
            N4 = (1.0 - xi) * (1.0 + eta) / 4.0;

            nodes[j * (m + 1) + i].ID = j * (m + 1) + i;

            nodes[j * (m + 1) + i].coord[0] = input.p1[0] * N1 + input.p2[0] * N2 + input.p3[0] * N3 + input.p4[0] * N4;
            nodes[j * (m + 1) + i].coord[1] = input.p1[1] * N1 + input.p2[1] * N2 + input.p3[1] * N3 + input.p4[1] * N4;
            nodes[j * (m + 1) + i].coord[2] =
                bathymetry_function(nodes[j * (m + 1) + i].coord[0], nodes[j * (m + 1) + i].coord[1]);
            ;
        }
    }

    std::vector<Element> elements(2 * m * n);

    switch (input.pattern) {
        case simple:
            simple_pattern_tri(m, n, elements);
            break;
        case zigzag:
            zigzag_pattern_tri(m, n, elements);
            break;
        case checker:
            checker_pattern_tri(m, n, elements);
            break;
    }

    std::vector<Boundary> boundaries(input.boundary);

    boundaries[0].nodes.resize(m + 1);
    boundaries[1].nodes.resize(n + 1);
    boundaries[2].nodes.resize(m + 1);
    boundaries[3].nodes.resize(n + 1);

    for (uint i = 0; i <= m; ++i) {
        boundaries[0].nodes[i] = i;
        boundaries[2].nodes[i] = i + (m + 1) * n;
    }

    for (uint i = 0; i <= n; ++i) {
        boundaries[1].nodes[i] = m + i * (m + 1);
        boundaries[3].nodes[i] = i * (m + 1);
    }

    std::ofstream file(std::string{input.mesh_name + ".14"});

    file << std::fixed << std::setprecision(12);
    file << "quadrilateral\n";
    file << 2 * m * n << "    " << (m + 1) * (n + 1) << '\n';

    for (auto& node : nodes) {
        file << node.ID << ' ';
        file << node.coord[0] << ' ';
        file << node.coord[1] << ' ';
        file << node.coord[2] << ' ';
        file << '\n';
    }

    for (auto& element : elements) {
        file << element.ID << ' ';
        file << element.type << ' ';
        file << element.nodes[0] << ' ';
        file << element.nodes[1] << ' ';
        file << element.nodes[2] << ' ';
        file << '\n';
    }

    uint n_land     = 0;
    uint n_tide     = 0;
    uint n_flow     = 0;
    uint n_function = 0;
    uint n_outflow  = 0;

    uint n_land_node     = 0;
    uint n_tide_node     = 0;
    uint n_flow_node     = 0;
    uint n_function_node = 0;
    uint n_outflow_node  = 0;

    for (uint n_bound = 0; n_bound < 4; ++n_bound) {
        if (boundaries[n_bound].type == 0) {
            n_land++;
            n_land_node += boundaries[n_bound].nodes.size();
        } else if (boundaries[n_bound].type == 1) {
            n_tide++;
            n_tide_node += boundaries[n_bound].nodes.size();
        } else if (boundaries[n_bound].type == 2) {
            n_flow++;
            n_flow_node += boundaries[n_bound].nodes.size();
        } else if (boundaries[n_bound].type == 3) {
            n_function++;
            n_function_node += boundaries[n_bound].nodes.size();
        } else if (boundaries[n_bound].type == 4) {
            n_outflow++;
            n_outflow_node += boundaries[n_bound].nodes.size();
        }
    }

    file << n_tide << " = Number of open boundaries\n";
    file << n_tide_node << " = Total number of open boundary nodes\n";

    uint i = 1;
    for (uint n_bound = 0; n_bound < 4; ++n_bound) {
        if (boundaries[n_bound].type == SWE::BoundaryTypes::tide) {
            file << boundaries[n_bound].nodes.size() << " = Number of nodes for open boundary " << i << '\n';

            std::for_each(boundaries[n_bound].nodes.begin(), boundaries[n_bound].nodes.end(), [&file](uint val) {
                file << val << '\n';
            });

            i++;
        }
    }

    file << n_land + n_flow + n_function + n_outflow << " = Number of land boundaries\n";
    file << n_land_node + n_flow_node + n_function_node + n_outflow_node << " = Total number of land boundary nodes\n";

    i = 1;
    for (uint n_bound = 0; n_bound < 4; ++n_bound) {
        if (boundaries[n_bound].type == SWE::BoundaryTypes::land) {
            file << boundaries[n_bound].nodes.size() << " 0 = Number of nodes for land boundary " << i << '\n';

            std::for_each(boundaries[n_bound].nodes.begin(), boundaries[n_bound].nodes.end(), [&file](uint val) {
                file << val << '\n';
            });

            i++;
        }
        if (boundaries[n_bound].type == SWE::BoundaryTypes::flow) {
            file << boundaries[n_bound].nodes.size() << " 2 = Number of nodes for land boundary " << i << '\n';

            std::for_each(boundaries[n_bound].nodes.begin(), boundaries[n_bound].nodes.end(), [&file](uint val) {
                file << val << '\n';
            });

            i++;
        }
        if (boundaries[n_bound].type == SWE::BoundaryTypes::function) {
            file << boundaries[n_bound].nodes.size() << " 77 = Number of nodes for land boundary " << i << '\n';

            std::for_each(boundaries[n_bound].nodes.begin(), boundaries[n_bound].nodes.end(), [&file](uint val) {
                file << val << '\n';
            });

            i++;
        }
        if (boundaries[n_bound].type == SWE::BoundaryTypes::outflow) {
            file << boundaries[n_bound].nodes.size() << " 88 = Number of nodes for land boundary " << i << '\n';

            std::for_each(boundaries[n_bound].nodes.begin(), boundaries[n_bound].nodes.end(), [&file](uint val) {
                file << val << '\n';
            });

            i++;
        }
    }

    file << "0 = Number of generic boundaries\n";
    file << "0 = Total number of generic boundary nodes\n";

    WriteBCISFile(input.mesh_name, boundaries);
}

void WriteBCISFile(const std::string& mesh_name, const std::vector<Boundary>& boundaries) {
    std::ofstream file(std::string{mesh_name + ".bcis"});

    for (const Boundary& boundary : boundaries) {
        if ((boundary.type == SWE::BoundaryTypes::tide) || (boundary.type == SWE::BoundaryTypes::flow)) {
            file << static_cast<int>(boundary.type) << ' ' << boundary.nodes.size() << '\n';
            // number of constituents currently hard-coded to 1
            file << 1 << '\n';
            file << boundary.frequency << ' ' << boundary.forcing_factor << ' ' << boundary.equilibrium_argument
                 << '\n';
            for (uint node_id : boundary.nodes) {
                file << node_id << ' ' << 1.0 << ' ' << 0.0 << '\n';
            }
        }
    }
}

MeshGeneratorInput::MeshGeneratorInput(const std::string& input_string) : boundary(4) {
    YAML::Node yaml_input = YAML::LoadFile(input_string);

    auto throw_missing_node = [](const std::string& str) {
        std::string err_msg{"Error: "};
        err_msg += str;
        err_msg += " yaml-node not specified\n";
        throw std::logic_error(err_msg);
    };

    if (yaml_input["x1"]) {
        this->p1[0] = yaml_input["x1"].as<double>();
    } else {
        throw_missing_node("x1");
    }

    if (yaml_input["y1"]) {
        this->p1[1] = yaml_input["y1"].as<double>();
    } else {
        throw_missing_node("y1");
    }

    if (yaml_input["x2"]) {
        this->p2[0] = yaml_input["x2"].as<double>();
    } else {
        throw_missing_node("x2");
    }

    if (yaml_input["y2"]) {
        this->p2[1] = yaml_input["y2"].as<double>();
    } else {
        throw_missing_node("y2");
    }

    if (yaml_input["x3"]) {
        this->p3[0] = yaml_input["x3"].as<double>();
    } else {
        throw_missing_node("x3");
    }

    if (yaml_input["y3"]) {
        this->p3[1] = yaml_input["y3"].as<double>();
    } else {
        throw_missing_node("y3");
    }

    if (yaml_input["x4"]) {
        this->p4[0] = yaml_input["x4"].as<double>();
    } else {
        throw_missing_node("x4");
    }

    if (yaml_input["y4"]) {
        this->p4[1] = yaml_input["y4"].as<double>();
    } else {
        throw_missing_node("y4");
    }

    if (yaml_input["num_x_subdivisions"]) {
        this->num_x_subdivisions = yaml_input["num_x_subdivisions"].as<uint>();
    } else {
        throw_missing_node("num_x_subdivisions");
    }

    if (yaml_input["num_y_subdivisions"]) {
        this->num_y_subdivisions = yaml_input["num_y_subdivisions"].as<uint>();
    } else {
        throw_missing_node("num_y_subdivisions");
    }

    if (yaml_input["pattern"]) {
        this->pattern = static_cast<Pattern>(yaml_input["pattern"].as<uint>());
        if (this->pattern > 2) {
            throw std::logic_error("Error: invalid pattern (" + std::to_string(this->pattern) + ") specified");
        }
    } else {
        throw_missing_node("pattern");
    }

    if (yaml_input["boundary"]) {
        YAML::Node yaml_boundaries = yaml_input["boundary"];
        if (!(yaml_boundaries.IsSequence() && yaml_boundaries.size() == 4)) {
            std::string err_msg{
                "Error: please check that the boundaries node is a "
                "sequence with 4 members"};
            throw std::logic_error(err_msg);
        }

        uint bid{0};
        for (auto it = yaml_boundaries.begin(); it != yaml_boundaries.end(); ++it) {
            Boundary& bd             = boundary[bid++];
            YAML::Node boundary_node = *it;

            std::string type_str = boundary_node["type"].as<std::string>();
            if (type_str == "land") {
                bd.type = SWE::BoundaryTypes::land;
            } else if (type_str == "tide") {
                bd.type = SWE::BoundaryTypes::tide;
            } else if (type_str == "flow") {
                bd.type = SWE::BoundaryTypes::flow;
            } else if (type_str == "function") {
                bd.type = SWE::BoundaryTypes::function;
            } else if (type_str == "outflow") {
                bd.type = SWE::BoundaryTypes::outflow;
            } else {
                std::string err_msg{"Error: Boundary type: " + type_str + " undefined"};
                throw std::logic_error(err_msg);
            }

            if (bd.type == SWE::BoundaryTypes::tide || bd.type == SWE::BoundaryTypes::flow) {
                if (!boundary_node["frequency"] || !boundary_node["forcing_factor"] ||
                    !boundary_node["equilibrium_argument"]) {
                    std::string err_msg{"Error: mal-formatted boundary node"};
                    throw std::logic_error(err_msg);
                }

                bd.frequency            = boundary_node["frequency"].as<double>();
                bd.forcing_factor       = boundary_node["forcing_factor"].as<double>();
                bd.equilibrium_argument = boundary_node["equilibrium_argument"].as<double>();
            }
        }

    } else {
        throw_missing_node("boundary");
    }

    this->mesh_name = yaml_input["mesh_name"] ? yaml_input["mesh_name"].as<std::string>() : "quad_mesh";

    this->Summarize();
}

void MeshGeneratorInput::Summarize() {
    std::cout << "MeshGenerator Input\n"
              << "  x1: " << this->p1[0] << '\n'
              << "  y1: " << this->p1[1] << '\n'
              << "  x2: " << this->p2[0] << '\n'
              << "  y2: " << this->p2[1] << '\n'
              << "  x3: " << this->p3[0] << '\n'
              << "  y3: " << this->p3[1] << '\n'
              << "  x4: " << this->p4[0] << '\n'
              << "  y4: " << this->p4[1] << '\n'
              << '\n'
              << "  num_x_subdivisions: " << num_x_subdivisions << '\n'
              << "  num_y_subdivisions: " << num_y_subdivisions << "\n\n";

    switch (this->pattern) {
        case simple:
            std::cout << "  pattern: simple\n\n";
            break;
        case zigzag:
            std::cout << "  pattern: zizag\n\n";
            break;
        case checker:
            std::cout << "  pattern: checker\n\n";
            break;
    }

    std::cout << "  boundary:\n";
    for (uint i = 0; i < this->boundary.size(); ++i) {
        std::cout << "  " << i << ":\n" << this->boundary[i];
    }
    std::cout << "\n\n";

    std::cout << "  mesh name: " << this->mesh_name << '\n';
}

void simple_pattern_tri(uint m, uint n, std::vector<Element>& elements) {
    // mesh with triangular elements simple pattern
    for (uint j = 0; j < n; ++j) {
        for (uint i = 0; i < m; ++i) {
            elements[2 * i + 2 * m * j].ID   = 2 * i + 2 * m * j;
            elements[2 * i + 2 * m * j].type = 3;
            elements[2 * i + 2 * m * j].nodes =
                std::vector<uint>{i + (m + 1) * j, i + m + 2 + (m + 1) * j, i + m + 1 + (m + 1) * j};

            elements[2 * i + 2 * m * j + 1].ID   = 2 * i + 2 * m * j + 1;
            elements[2 * i + 2 * m * j + 1].type = 3;
            elements[2 * i + 2 * m * j + 1].nodes =
                std::vector<uint>{i + (m + 1) * j, i + 1 + (m + 1) * j, i + m + 2 + (m + 1) * j};
        }
    }
}

void zigzag_pattern_tri(uint m, uint n, std::vector<Element>& elements) {
    // mesh with triangular elements zigzag pattern
    for (uint j = 0; j < n; j += 2) {
        for (uint i = 0; i < m; ++i) {
            elements[2 * i + 2 * m * j].ID   = 2 * i + 2 * m * j;
            elements[2 * i + 2 * m * j].type = 3;
            elements[2 * i + 2 * m * j].nodes =
                std::vector<uint>{i + (m + 1) * j, i + m + 2 + (m + 1) * j, i + m + 1 + (m + 1) * j};

            elements[2 * i + 2 * m * j + 1].ID   = 2 * i + 2 * m * j + 1;
            elements[2 * i + 2 * m * j + 1].type = 3;
            elements[2 * i + 2 * m * j + 1].nodes =
                std::vector<uint>{i + (m + 1) * j, i + 1 + (m + 1) * j, i + m + 2 + (m + 1) * j};
        }
    }

    for (uint j = 1; j < n; j += 2) {
        for (uint i = 0; i < m; ++i) {
            elements[2 * i + 2 * m * j].ID   = 2 * i + 2 * m * j;
            elements[2 * i + 2 * m * j].type = 3;
            elements[2 * i + 2 * m * j].nodes =
                std::vector<uint>{i + (m + 1) * j, i + 1 + (m + 1) * j, i + m + 1 + (m + 1) * j};

            elements[2 * i + 2 * m * j + 1].ID   = 2 * i + 2 * m * j + 1;
            elements[2 * i + 2 * m * j + 1].type = 3;
            elements[2 * i + 2 * m * j + 1].nodes =
                std::vector<uint>{i + 1 + (m + 1) * j, i + m + 2 + (m + 1) * j, i + m + 1 + (m + 1) * j};
        }
    }
}

void checker_pattern_tri(uint m, uint n, std::vector<Element>& elements) {
    // mesh with triangular elements checker pattern
    for (uint j = 0; j < n; ++j) {
        for (uint i = j % 2; i < m; i += 2) {
            elements[2 * i + 2 * m * j].ID   = 2 * i + 2 * m * j;
            elements[2 * i + 2 * m * j].type = 3;
            elements[2 * i + 2 * m * j].nodes =
                std::vector<uint>{i + (m + 1) * j, i + m + 2 + (m + 1) * j, i + m + 1 + (m + 1) * j};

            elements[2 * i + 2 * m * j + 1].ID   = 2 * i + 2 * m * j + 1;
            elements[2 * i + 2 * m * j + 1].type = 3;
            elements[2 * i + 2 * m * j + 1].nodes =
                std::vector<uint>{i + (m + 1) * j, i + 1 + (m + 1) * j, i + m + 2 + (m + 1) * j};
        }
    }

    for (uint j = 0; j < n; ++j) {
        for (uint i = (j + 1) % 2; i < m; i += 2) {
            elements[2 * i + 2 * m * j].ID   = 2 * i + 2 * m * j;
            elements[2 * i + 2 * m * j].type = 3;
            elements[2 * i + 2 * m * j].nodes =
                std::vector<uint>{i + (m + 1) * j, i + 1 + (m + 1) * j, i + m + 1 + (m + 1) * j};

            elements[2 * i + 2 * m * j + 1].ID   = 2 * i + 2 * m * j + 1;
            elements[2 * i + 2 * m * j + 1].type = 3;
            elements[2 * i + 2 * m * j + 1].nodes =
                std::vector<uint>{i + 1 + (m + 1) * j, i + m + 2 + (m + 1) * j, i + m + 1 + (m + 1) * j};
        }
    }
}
