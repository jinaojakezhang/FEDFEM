/*
MIT License

Copyright (c) 2021 Jinao Zhang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <omp.h>
using namespace    std;
static const int   NUM_THREADS(omp_get_max_threads());
static const float ABS_ZERO_T(-273.15f),
                   STEFAN_BOLTZMANN(5.67f*0.00000001f);

// matrix computation/operations (mat: matrix, 33: 3 rows by 3 columns, x: multiplication, T: transpose, Det: determinant, Inv: inverse)
void mat33x34    (const float A[3][3], const float B[3][4], float AB[3][4]);
void mat34x34T   (const float A[3][4], const float B[3][4], float AB[3][3]);
void mat34Tx34   (const float A[3][4], const float B[3][4], float AB[4][4]);
void mat44xScalar(const float A[4][4], const float b,       float Ab[4][4]);
void matDet33    (const float A[3][3], float &detA);
void matInv33    (const float A[3][3], float invA[3][3], float &detA);

// classes
class Node;
class T4;
class Model;
class ModelStates;

// methods
Model*       readModel       (int argc, char **argv);
void         printInfo       (const Model& model);
ModelStates* runSimulation   (const Model& model);
void         initBC          (const Model& model, ModelStates& modelstates);
void         computeRunTimeBC(const Model& model, ModelStates& modelstates, const size_t curr_step);
bool         computeOneStep  (const Model& model, ModelStates& modelstates);
int          exportVTK       (const Model& model, const ModelStates& modelstates);

class Node
{
public:
    const unsigned int m_idx;
    const float        m_x, m_y, m_z;
    Node(const unsigned int idx, const float x, const float y, const float z) :
        m_idx(idx), m_x(x), m_y(y), m_z(z) {};
};

class T4
{
public:
    const unsigned int m_idx, m_n_idx[4];
    const float        m_DHDr[3][4];
    float              m_DHDX[3][4],
                       m_D[3][3], m_K[4][4], // D: conductivity; K: conduction
                       m_Vol, m_mass;
    const string       m_material_type;
    vector<float>      m_material_vals;
    T4(const unsigned int idx, const Node& n1, const Node& n2, const Node& n3, const Node& n4, const float rho, const string material_type, const vector<float>& material_vals) :
        m_idx(idx), m_n_idx{ n1.m_idx, n2.m_idx, n3.m_idx, n4.m_idx }, m_material_type(material_type), m_material_vals(material_vals.cbegin(), material_vals.cend()),
        m_DHDr{{-1, 1, 0, 0},
               {-1, 0, 1, 0},
               {-1, 0, 0, 1}}
    {
        float n_coords[3][4], J0[3][3], detJ0(0.f), invJ0[3][3];
        n_coords[0][0] = n1.m_x; n_coords[1][0] = n1.m_y; n_coords[2][0] = n1.m_z;
        n_coords[0][1] = n2.m_x; n_coords[1][1] = n2.m_y; n_coords[2][1] = n2.m_z;
        n_coords[0][2] = n3.m_x; n_coords[1][2] = n3.m_y; n_coords[2][2] = n3.m_z;
        n_coords[0][3] = n4.m_x; n_coords[1][3] = n4.m_y; n_coords[2][3] = n4.m_z;
        mat34x34T(m_DHDr, n_coords, J0);
        matInv33(J0, invJ0, detJ0);
        m_Vol = detJ0 / 6.f;
        m_mass = rho * m_Vol;
        mat33x34(invJ0, m_DHDr, m_DHDX);
        if (m_material_type == "T_ISO") // [0]=c, [1]=k
        {
            memset(m_D, 0, sizeof(float) * 3 * 3);
            m_D[0][0] = m_material_vals[1]; m_D[1][1] = m_material_vals[1]; m_D[2][2] = m_material_vals[1];
        }
        else if (m_material_type == "T_ORTHO") // [0]=c, [1]=k11, [2]=k22, [3]=k33
        {
            memset(m_D, 0, sizeof(float) * 3 * 3);
            m_D[0][0] = m_material_vals[1]; m_D[1][1] = m_material_vals[2]; m_D[2][2] = m_material_vals[3];
        }
        else if (m_material_type == "T_ANISO") // [0]=c, [1]=k11, [2]=k12, [3]=k13, [4]=k22, [5]=k23, [6]=k33
        {
            memset(m_D, 0, sizeof(float) * 3 * 3);
            m_D[0][0] = m_material_vals[1]; m_D[0][1] = m_material_vals[2]; m_D[0][2] = m_material_vals[3];
            m_D[1][0] = m_D[0][1];          m_D[1][1] = m_material_vals[4]; m_D[1][2] = m_material_vals[5];
            m_D[2][0] = m_D[0][2];          m_D[2][1] = m_D[1][2];          m_D[2][2] = m_material_vals[6];
        }
        float temp[3][4];
        mat33x34(m_D, m_DHDX, temp);
        mat34Tx34(m_DHDX, temp, m_K);
        mat44xScalar(m_K, m_Vol, m_K);
    };
};

class Model
{
public:
    vector<Node*>        m_nodes;
    vector<T4*>          m_tets;
    size_t               m_num_BCs,    m_num_steps,  m_num_DOFs;
    vector<unsigned int> m_hflux_idx,  m_convc_idx,  m_radia_idx, m_fixT_idx, m_bhflux_idx;
    vector<float>        m_hflux_mag,
                         m_convc_coef, m_convc_surf, m_convc_refT, m_convc_const1,
                         m_radia_emis, m_radia_surf, m_radia_refT, m_radia_const1, m_radia_const2,
                         m_fixT_mag,
                         m_bhflux_mag,
                         m_material_vals;
    float                m_dt, m_total_t, m_T0, m_rho;
    const string         m_fname;
    string               m_ele_type, m_material_type;
    unsigned int         m_node_begin_index, m_ele_begin_index,
                        *m_ele_node_local_idx_pair,
                        *m_tracking_num_eles_i_eles_per_node_j;
    Model(const string fname) :
        m_nodes     (0), m_tets       (0),
        m_num_BCs   (0), m_num_steps  (0), m_num_DOFs  (0),
        m_hflux_idx (0), m_hflux_mag  (0),
        m_convc_idx (0), m_convc_coef (0), m_convc_surf(0), m_convc_refT(0), m_convc_const1(0),
        m_radia_idx (0), m_radia_emis (0), m_radia_surf(0), m_radia_refT(0), m_radia_const1(0), m_radia_const2(0),
        m_fixT_idx  (0), m_fixT_mag   (0),
        m_bhflux_idx(0), m_bhflux_mag (0),
        m_material_vals(0), m_dt(0.f), m_total_t(0.f), m_T0(0.f), m_rho(0.f),
        m_fname(fname), m_ele_type(""), m_material_type(""),
        m_node_begin_index(0), m_ele_begin_index(0),
        m_ele_node_local_idx_pair(nullptr), m_tracking_num_eles_i_eles_per_node_j(nullptr) {};
    ~Model()
    {
        for (Node* node : m_nodes) { delete node; }
        for (T4* tet : m_tets)     { delete tet; }
        delete[] m_ele_node_local_idx_pair;
        delete[] m_tracking_num_eles_i_eles_per_node_j;
    };
    void postCreate()
    {
        // below: provide indexing for nodal states (e.g., individual ele nodal thermal loads) to avoid race condition in parallel computing
        m_tracking_num_eles_i_eles_per_node_j = new unsigned int[m_nodes.size() * 2]; memset(m_tracking_num_eles_i_eles_per_node_j, 0, sizeof(unsigned int) * m_nodes.size() * 2);
        vector<vector<unsigned int>> nodes_ele_node_local_idx_pair(m_nodes.size());
        for (T4* tet : m_tets) { for (unsigned int m = 0; m < 4; m++) { nodes_ele_node_local_idx_pair[tet->m_n_idx[m]].push_back(tet->m_idx); nodes_ele_node_local_idx_pair[tet->m_n_idx[m]].push_back(m); } }
        vector<unsigned int> eles_per_node(m_nodes.size(), 0);
        unsigned int length(0);
        for (size_t m = 0; m < m_nodes.size(); m++)
        {
            eles_per_node[m] = (unsigned int)nodes_ele_node_local_idx_pair[m].size() / 2;
            length += eles_per_node[m];
        }
        m_ele_node_local_idx_pair = new unsigned int[length * 2];
        unsigned int* p_ele_node_local_idx_pair = m_ele_node_local_idx_pair, tracking(0);
        for (size_t m = 0; m < m_nodes.size(); m++)
        {
            m_tracking_num_eles_i_eles_per_node_j[m * 2 + 0] = tracking;
            m_tracking_num_eles_i_eles_per_node_j[m * 2 + 1] = eles_per_node[m]; tracking += eles_per_node[m];
            for (size_t n = 0; n < eles_per_node[m]; n++)
            {
                *p_ele_node_local_idx_pair = nodes_ele_node_local_idx_pair[m][n * 2 + 0]; p_ele_node_local_idx_pair++; // tet->m_idx
                *p_ele_node_local_idx_pair = nodes_ele_node_local_idx_pair[m][n * 2 + 1]; p_ele_node_local_idx_pair++; // m
            }
        }
    }
};

class ModelStates
{
public:
    vector<float> m_external_Q, m_external_Q0, m_ele_nodal_internal_Q, // individual ele nodal internal Q to avoid race condition, can be summed to get internal_Q for nodes
                  m_fixT_mag,
                  m_constA,
                  m_curr_T,     m_next_T;
    vector<bool>  m_fixT_flag;
    ModelStates(const Model& model) :
        m_external_Q(model.m_num_DOFs,        0.f), m_external_Q0(model.m_num_DOFs,        0.f), m_ele_nodal_internal_Q(model.m_tets.size() * 4, 0.f), // each tet has 4 nodes, each node has 1 (T) DOF
        m_fixT_mag  (model.m_num_DOFs,        0.f),
        m_constA    (model.m_num_DOFs,        0.f),
        m_curr_T    (model.m_num_DOFs, model.m_T0), m_next_T     (model.m_num_DOFs, model.m_T0),
        m_fixT_flag (model.m_num_DOFs,      false)
    {
        vector<float> nodal_mass(model.m_num_DOFs, 0.f);
        for (T4* tet : model.m_tets) { for (size_t m = 0; m < 4; m++) { nodal_mass[tet->m_n_idx[m]] += tet->m_mass / 4.f; } }
        for (size_t i = 0; i < model.m_num_DOFs; i++) { m_constA[i] = model.m_dt / (nodal_mass[i] * model.m_material_vals[0]); }
    };
};

int main(int argc, char **argv)
{
    Model* model = readModel(argc, argv);
    if (model != nullptr)
    {
        printInfo(*model);
        ModelStates* modelstates = runSimulation(*model);
        if (modelstates != nullptr)
        {
            int exit = exportVTK(*model, *modelstates);
            delete model;
            delete modelstates;
            return exit;
        }
        else { return EXIT_FAILURE; }
    }
    else { return EXIT_FAILURE; }
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
Model* readModel(int argc, char **argv)
{
    if (argc - 1 == 0) { cerr << "\n\tError: missing input argument (e.g., T_Iso.txt)." << endl; return nullptr; }
    FILE* file;
    if (fopen_s(&file, argv[1], "r") != 0) { cerr << "\n\tError: cannot open file: " << argv[1] << endl; return nullptr; }
    else
    {
        Model* model = new Model(argv[1]);
        char buffer[256];
        unsigned int idx(0); float x(0.f), y(0.f), z(0.f);
        fscanf_s(file, "%u %f %f %f", &idx, &x, &y, &z); model->m_node_begin_index = idx; model->m_nodes.push_back(new Node(idx - model->m_node_begin_index, x, y, z)); // for first node only, to get node begin index
        while (fscanf_s(file, "%u %f %f %f", &idx, &x, &y, &z)) { model->m_nodes.push_back(new Node(idx - model->m_node_begin_index, x, y, z)); } // internally, node index starts at 0
        fscanf_s(file, "%s", buffer, (unsigned int)sizeof(buffer)); model->m_material_type = buffer;
        if      (model->m_material_type == "T_ISO")   { float c(0.f), k(0.f); fscanf_s(file, "%f %f", &c, &k); model->m_material_vals.push_back(c); model->m_material_vals.push_back(k); }
        else if (model->m_material_type == "T_ORTHO") { float c(0.f), k11(0.f), k22(0.f), k33(0.f); fscanf_s(file, "%f %f %f %f", &c, &k11, &k22, &k33); model->m_material_vals.push_back(c); model->m_material_vals.push_back(k11); model->m_material_vals.push_back(k22); model->m_material_vals.push_back(k33); }
        else if (model->m_material_type == "T_ANISO") { float c(0.f), k11(0.f), k12(0.f), k13(0.f), k22(0.f), k23(0.f), k33(0.f); fscanf_s(file, "%f %f %f %f %f %f %f", &c, &k11, &k12, &k13, &k22, &k23, &k33); model->m_material_vals.push_back(c); model->m_material_vals.push_back(k11); model->m_material_vals.push_back(k12); model->m_material_vals.push_back(k13); model->m_material_vals.push_back(k22); model->m_material_vals.push_back(k23); model->m_material_vals.push_back(k33); }
        else if (model->m_material_type == "other_conductivity_types") { /*add your code here*/ }
        fscanf_s(file, "%s %f", buffer, (unsigned int)sizeof(buffer), &model->m_rho);
        fscanf_s(file, "%s", buffer, (unsigned int)sizeof(buffer)); model->m_ele_type = buffer;
        unsigned int n1_idx(0), n2_idx(0), n3_idx(0), n4_idx(0);
        fscanf_s(file, "%u %u %u %u %u", &idx, &n1_idx, &n2_idx, &n3_idx, &n4_idx); model->m_ele_begin_index = idx; model->m_tets.push_back(new T4(idx - model->m_ele_begin_index,
                                                                                                                                                   *model->m_nodes[n1_idx - model->m_node_begin_index], *model->m_nodes[n2_idx - model->m_node_begin_index], *model->m_nodes[n3_idx - model->m_node_begin_index], *model->m_nodes[n4_idx - model->m_node_begin_index],
                                                                                                                                                    model->m_rho, model->m_material_type, model->m_material_vals)); // for first ele only, to get ele begin index
        while (fscanf_s(file, "%u %u %u %u %u", &idx, &n1_idx, &n2_idx, &n3_idx, &n4_idx)) { model->m_tets.push_back(new T4(idx - model->m_ele_begin_index,
                                                                                                                            *model->m_nodes[n1_idx - model->m_node_begin_index], *model->m_nodes[n2_idx - model->m_node_begin_index], *model->m_nodes[n3_idx - model->m_node_begin_index], *model->m_nodes[n4_idx - model->m_node_begin_index],
                                                                                                                             model->m_rho, model->m_material_type, model->m_material_vals)); } // internally, ele index starts at 0
        while (fscanf_s(file, "%s", buffer, (unsigned int)sizeof(buffer)))
        {
            string BC_type(buffer);
            if (BC_type == "<HFlux>") // Nodal heat flux
            {
                float q(0.f); fscanf_s(file, "%f", &q);
                while (fscanf_s(file, "%u", &idx)) { model->m_hflux_idx.push_back(idx - model->m_node_begin_index); model->m_hflux_mag.push_back(q); }
                model->m_num_BCs++;
            }
            else if (BC_type == "<Convc>") // Convection
            {
                float coef(0.f), surf(0.f), refT(0.f); fscanf_s(file, "%f %f %f", &coef, &surf, &refT);
                while (fscanf_s(file, "%u", &idx)) { model->m_convc_idx.push_back(idx - model->m_node_begin_index); model->m_convc_coef.push_back(coef); model->m_convc_surf.push_back(surf); model->m_convc_refT.push_back(refT); model->m_convc_const1.push_back(coef * surf); }
                model->m_num_BCs++;
            }
            else if (BC_type == "<Radia>") // Radiation
            {
                float emis(0.f), surf(0.f), refT(0.f); fscanf_s(file, "%f %f %f", &emis, &surf, &refT);
                while (fscanf_s(file, "%u", &idx)) { model->m_radia_idx.push_back(idx - model->m_node_begin_index); model->m_radia_emis.push_back(emis); model->m_radia_surf.push_back(surf); model->m_radia_refT.push_back(refT); model->m_radia_const1.push_back(STEFAN_BOLTZMANN * emis * surf); model->m_radia_const2.push_back(powf(refT - ABS_ZERO_T, 4.f)); }
                model->m_num_BCs++;
            }
            else if (BC_type == "<FixT>") // Fixed temperature
            {
                float constT(0.f); fscanf_s(file, "%f", &constT);
                while (fscanf_s(file, "%u", &idx)) { model->m_fixT_idx.push_back(idx - model->m_node_begin_index); model->m_fixT_mag.push_back(constT); }
                model->m_num_BCs++;
            }
            else if (BC_type == "<BodyHFlux>") // Body heat flux
            {
                float q(0.f); fscanf_s(file, "%f", &q);
                vector<float> nodal_q(model->m_nodes.size(), 0.f);
                while (fscanf_s(file, "%u", &idx)) { for (size_t m = 0; m < 4; m++) { nodal_q[model->m_tets[idx - model->m_ele_begin_index]->m_n_idx[m]] += q * model->m_tets[idx - model->m_ele_begin_index]->m_Vol / 4.f; } }
                for (unsigned int i = 0; i < nodal_q.size(); i++) { if (nodal_q[i] != 0) { model->m_bhflux_idx.push_back(i); model->m_bhflux_mag.push_back(nodal_q[i]); } }
                model->m_num_BCs++;
            }
            else if (BC_type == "other_BC_types") { /*add your code here*/ }
            else if (BC_type == "</BC>") { break; }
        }
        fscanf_s(file, "%s %f", buffer, (unsigned int)sizeof(buffer), &model->m_T0);
        fscanf_s(file, "%s %f", buffer, (unsigned int)sizeof(buffer), &model->m_dt);
        fscanf_s(file, "%s %f", buffer, (unsigned int)sizeof(buffer), &model->m_total_t);
        fclose(file);
        model->m_num_steps = (size_t)ceil(model->m_total_t / model->m_dt);
        model->m_num_DOFs = model->m_nodes.size() * 1;
        model->postCreate();
        return model;
    }
}

void printInfo(const Model& model)
{
    cout << endl;
    cout << "\t---------------------------------------------------------------------------------------------------" << endl;
    cout << "\t| Oepn-source (OpenMP) implmentation of:                                                          |" << endl;
    cout << "\t|                   <Fast explicit dynamics finite element algorithm for transient heat transfer. |" << endl;
    cout << "\t|                                                               Zhang, J., & Chauhan, S. (2019).  |" << endl;
    cout << "\t|                                       International Journal of Thermal Sciences, 139, 160-175.  |" << endl;
    cout << "\t|                                                         doi:10.1016/j.ijthermalsci.2019.01.030> |" << endl;
    cout << "\t|                                                                                  by Jinao Zhang |" << endl;
    cout << "\t---------------------------------------------------------------------------------------------------" << endl;
    cout << "\tModel:\t\t"      << model.m_fname.c_str()         << endl;
    cout << "\tNodes:\t\t"      << model.m_nodes.size()          << " (" << model.m_num_DOFs    << " DOFs)" << endl;
    cout << "\tElements:\t"     << model.m_tets.size()           << " (" << model.m_ele_type.c_str() << ")" << endl;
    cout << "\tEleMaterial:\t"  << model.m_material_type.c_str() << ":"; for (const float val : model.m_material_vals) { cout << " " << val; } cout << "; Density: " << model.m_rho << endl;
    cout << "\tBC:\t\t"         << model.m_num_BCs               << endl;
    cout << "\tInitialTemp.:\t" << model.m_T0                    << endl;
    cout << "\tTimeStep:\t"     << model.m_dt                    << endl;
    cout << "\tTotalTime:\t"    << model.m_total_t               << endl;
    cout << "\tNumSteps:\t"     << model.m_num_steps             << endl;
    cout << "\n\tNode index starts at " << model.m_node_begin_index << "." << endl;
    cout << "  \tElem index starts at " << model.m_ele_begin_index  << "." << endl;
}

ModelStates* runSimulation(const Model& model)
{
    ModelStates* modelstates = new ModelStates(model);
    initBC(model, *modelstates);
    size_t progress(0);
    auto start_t = chrono::high_resolution_clock::now();
    cout << "\n\tusing " << NUM_THREADS << " threads" << endl;
    cout << "\tcomputing..." << endl;
    for (size_t step = 0; step < model.m_num_steps; step++) // simulation loop
    {
        if ((float)(step + 1) / (float)model.m_num_steps * 100.f >= progress + 10) { progress += 10; cout << "\t\t\t(" << progress << "%)" << endl; }
        computeRunTimeBC(model, *modelstates, step);
        bool no_err = computeOneStep(model, *modelstates);
        if (!no_err) { delete modelstates; return nullptr; }
    }
    auto elapsed = chrono::high_resolution_clock::now() - start_t;
    long long t = chrono::duration_cast<chrono::milliseconds>(elapsed).count();
    cout << "\n\tComputation time:\t" << t << " ms" << endl;
    return modelstates;
}

void initBC(const Model& model, ModelStates& modelstates)
{
    fill(modelstates.m_external_Q.begin(),  modelstates.m_external_Q.end(),  0.f);
    fill(modelstates.m_external_Q0.begin(), modelstates.m_external_Q0.end(), 0.f);
    // BC:HFlux
    for (size_t i = 0; i < model.m_hflux_idx.size(); i++) { modelstates.m_external_Q0[model.m_hflux_idx[i]] += model.m_hflux_mag[i]; }
    // BC:BodyHFlux
    for (size_t i = 0; i < model.m_bhflux_idx.size(); i++) { modelstates.m_external_Q0[model.m_bhflux_idx[i]] += model.m_bhflux_mag[i]; }
    // BC:FixT
    for (size_t i = 0; i < model.m_fixT_idx.size(); i++) { modelstates.m_fixT_flag[model.m_fixT_idx[i]] = true; modelstates.m_fixT_mag[model.m_fixT_idx[i]] = model.m_fixT_mag[i]; }
    modelstates.m_external_Q = modelstates.m_external_Q0;
}

void computeRunTimeBC(const Model& model, ModelStates& modelstates, const size_t curr_step)
{
    // BC:Convc
    for (size_t i = 0; i < model.m_convc_idx.size(); i++) { modelstates.m_external_Q[model.m_convc_idx[i]] = modelstates.m_external_Q0[model.m_convc_idx[i]] - model.m_convc_const1[i] * (modelstates.m_curr_T[model.m_convc_idx[i]] - model.m_convc_refT[i]); }
    // BC:Radia
    for (size_t i = 0; i < model.m_radia_idx.size(); i++)
    {
        const float pow1 = powf(modelstates.m_curr_T[model.m_radia_idx[i]] - ABS_ZERO_T, 4.f);
        modelstates.m_external_Q[model.m_radia_idx[i]] = modelstates.m_external_Q0[model.m_radia_idx[i]] - model.m_radia_const1[i] * (pow1 - model.m_radia_const2[i]);
    }
}

bool computeOneStep(const Model& model, ModelStates& modelstates)
{
    bool no_err(true);
#pragma omp parallel num_threads(NUM_THREADS)
    {
        int id = omp_get_thread_num();
        for (int i = id; i < model.m_tets.size(); i += NUM_THREADS) // loop through tets to compute for thermal load contributions
        {
            T4 *tet = model.m_tets[i];
            // compute ele q
            for (size_t m = 0; m < 4; m++) { modelstates.m_ele_nodal_internal_Q[tet->m_idx * 4 + m] = tet->m_K[m][0] * modelstates.m_curr_T[tet->m_n_idx[0]]
                                                                                                    + tet->m_K[m][1] * modelstates.m_curr_T[tet->m_n_idx[1]]
                                                                                                    + tet->m_K[m][2] * modelstates.m_curr_T[tet->m_n_idx[2]]
                                                                                                    + tet->m_K[m][3] * modelstates.m_curr_T[tet->m_n_idx[3]]; }
        }
#pragma omp barrier
        for (int i = id; i < model.m_nodes.size(); i += NUM_THREADS) // loop through nodes to compute for new temperatures T
        {
            // assemble nodal thermal loads from individual ele nodal thermal loads, due to avoiding race condition
            float nodal_internal_Q(0.f);
            unsigned int tracking_num_eles(model.m_tracking_num_eles_i_eles_per_node_j[i * 2 + 0]),
                         eles_per_node    (model.m_tracking_num_eles_i_eles_per_node_j[i * 2 + 1]),
                         ele_idx(0), node_local_idx(0);
            for (unsigned int j = 0; j < eles_per_node; j++)
            {
                ele_idx        = model.m_ele_node_local_idx_pair[(tracking_num_eles + j) * 2 + 0];
                node_local_idx = model.m_ele_node_local_idx_pair[(tracking_num_eles + j) * 2 + 1];
                nodal_internal_Q += modelstates.m_ele_nodal_internal_Q[ele_idx * 4 + node_local_idx];
            }
            if (modelstates.m_fixT_flag[i] == true) { modelstates.m_next_T[i] = modelstates.m_fixT_mag[i]; } // apply BC:FixT
            else                                                                                             // explicit time integration
            {
                modelstates.m_next_T[i] = modelstates.m_curr_T[i] + modelstates.m_constA[i] * (modelstates.m_external_Q[i] - nodal_internal_Q);
                if (isnan(modelstates.m_next_T[i])) { no_err = false; }
            }
        }
    }
    if (!no_err) { cerr << "\n\tError: solution diverged, simulation aborted. Try a smaller time step." << endl; }
    else { modelstates.m_curr_T.swap(modelstates.m_next_T); }
    return no_err;
}

int exportVTK(const Model& model, const ModelStates& modelstates)
{
    const vector<string> outputs{ "T.vtk" }; // other outputs can be added by the user
    cout << "\n\texporting..." << endl;
    for (string vtk : outputs)
    {
        ofstream fout(vtk.c_str());
        if (fout.is_open())
        {
            fout << "# vtk DataFile Version 3.8" << endl;
            fout << vtk.c_str() << endl;
            fout << "ASCII" << endl;
            fout << "DATASET UNSTRUCTURED_GRID" << endl;
            fout << "POINTS " << model.m_nodes.size() << " float" << endl;
            for (Node* node : model.m_nodes) { fout << node->m_x << " " << node->m_y << " " << node->m_z << endl; }
            fout << "CELLS " << model.m_tets.size() << " " << model.m_tets.size() * (4 + 1) << endl;
            for (T4* tet : model.m_tets) { fout << 4 << " " << tet->m_n_idx[0] << " " << tet->m_n_idx[1] << " " << tet->m_n_idx[2] << " " << tet->m_n_idx[3] << endl; }
            fout << "CELL_TYPES " << model.m_tets.size() << endl;
            for (size_t i = 0; i < model.m_tets.size(); i++) { fout << 10 << endl; }
            fout << "POINT_DATA " << model.m_nodes.size() << endl;
            fout << "SCALARS " << vtk.c_str() << " float" << endl;
            fout << "LOOKUP_TABLE default" << endl;
            for (Node* node : model.m_nodes) { fout << modelstates.m_curr_T[node->m_idx] << endl; }
            cout << "\t\t\t" << vtk.c_str() << endl;
        }
        else { cerr << "\n\tError: cannot open " << vtk.c_str() << " for writing, results not saved." << endl; return EXIT_FAILURE; }
    }
    cout << "\tVTK saved." << endl;
    return EXIT_SUCCESS;
}

void mat33x34(const float A[3][3], const float B[3][4], float AB[3][4])
{
    memset(AB, 0, sizeof(float) * 3 * 4);
    for (size_t i = 0; i < 3; i++) { for (size_t j = 0; j < 4; j++) { for (size_t k = 0; k < 3; k++) { AB[i][j] += A[i][k] * B[k][j]; } } }
}
void mat34x34T(const float A[3][4], const float B[3][4], float AB[3][3])
{
    memset(AB, 0, sizeof(float) * 3 * 3);
    for (size_t i = 0; i < 3; i++) { for (size_t j = 0; j < 3; j++) { for (size_t k = 0; k < 4; k++) { AB[i][j] += A[i][k] * B[j][k]; } } }
}
void mat34Tx34(const float A[3][4], const float B[3][4], float AB[4][4])
{
    memset(AB, 0, sizeof(float) * 4 * 4);
    for (size_t i = 0; i < 4; i++) { for (size_t j = 0; j < 4; j++) { for (size_t k = 0; k < 3; k++) { AB[i][j] += A[k][i] * B[k][j]; } } }
}
void mat44xScalar(const float A[4][4], const float b, float Ab[4][4])
{
    for (size_t i = 0; i < 4; i++) { for (size_t j = 0; j < 4; j++) { Ab[i][j] = A[i][j] * b; } }
}
void matDet33(const float A[3][3], float &detA)
{
    detA = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[1][0] * (A[0][1] * A[2][2] - A[0][2] * A[2][1]) + A[2][0] * (A[0][1] * A[1][2] - A[0][2] * A[1][1]);
}
void matInv33(const float A[3][3], float invA[3][3], float &detA)
{
    matDet33(A, detA);
    invA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / detA; invA[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / detA; invA[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / detA;
    invA[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / detA; invA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / detA; invA[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / detA;
    invA[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / detA; invA[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / detA; invA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / detA;
}