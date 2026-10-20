#include "GMSTClustering.hh"
#include "Repulsion.hh"
#include "AAMCConstants.hh"
#include <algorithm>
#include <cmath>
#include <stack>

GMSTCluster::GMSTCluster(G4int Z_in, G4int A_in) : A(A_in), Z(Z_in){};

GMSTCluster::~GMSTCluster(){};

GMSTClustering::GMSTClustering() { CritDist = 100; };

GMSTClustering::GMSTClustering(G4double CD_in, G4double Aa_in, G4double Ab_in) {
  CritDist = CD_in;
  SpecAa = Aa_in;
  SpecAb = Ab_in;
};

GMSTClustering::~GMSTClustering(){};

Graph GMSTClustering::ClusterToGraph(aamcc::NucleonVector *nucleons, double V0pass) {
  G4int N = (G4int)nucleons->size();
  Graph g(N, (N > 1) ? N * (N - 1) / 2 : 0);

  // Узлы-нуклоны (A=1, масса = nucleonAverMass)
  std::vector<Node> nodes(N);
  for (G4int i = 0; i < N; i++) {
    aamcc::Nucleon *n = &(nucleons->at(i));
    nodes[i] = { n->GetX(), n->GetY(), n->GetZ(),
                 n->GetPx(), n->GetPy(), n->GetPz(),
                 nucleonAverMass, 1 };
  }

  for (G4int i = 0; i < N; i++) {
    for (G4int j = i + 1; j < N; j++) {
      g.addEdge(i, j, PairMetric({nodes[i]}, {nodes[j]}, this->CritDist, V0pass));
    }
  }
  return g;
}

double GMSTClustering::PairMetric(const std::vector<Node>& a, const std::vector<Node>& b,
                                  double dmin, double V0pass) const {
  // Центроиды групп
  double cxA = 0, cyA = 0, czA = 0;
  for (const auto& n : a) { cxA += n.x; cyA += n.y; czA += n.z; }
  const double invA = 1.0 / static_cast<double>(a.size());
  cxA *= invA; cyA *= invA; czA *= invA;

  double cxB = 0, cyB = 0, czB = 0;
  for (const auto& n : b) { cxB += n.x; cyB += n.y; czB += n.z; }
  const double invB = 1.0 / static_cast<double>(b.size());
  cxB *= invB; cyB *= invB; czB *= invB;

  const double dx = cxA - cxB;
  const double dy = cyA - cyB;
  const double dz = czA - czB;

  // Режим 0 — чисто пространственный
  if (metricMode == 0) {
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  // Если порог скорости выключен — ведём себя как режим 0
  if (V0pass <= 0.0) {
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  // Полный 4-импульс объединения a+b
  double Ptx = 0, Pty = 0, Ptz = 0, Etot = 0;
  for (const auto& n : a) {
    Ptx += n.px; Pty += n.py; Ptz += n.pz;
    Etot += std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
  }
  for (const auto& n : b) {
    Ptx += n.px; Pty += n.py; Ptz += n.pz;
    Etot += std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
  }

  // Буст в СЦМ объединения (вдоль z)
  double beta = Ptz / Etot;
  if (beta >=  1.0) beta =  1.0 - 1e-12;
  if (beta <= -1.0) beta = -1.0 + 1e-12;
  const double g = 1.0 / std::sqrt(1.0 - beta * beta);

  // DCM-подобный per-particle cut: максимальная скорость составляющей в СЦМ объединения
  double d_V_max = 0;
  for (const auto& n : a) {
    const double e = std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
    const double pz_p = g * (n.pz - beta * e);
    const double e_p  = g * (e - beta * n.pz);
    const double vx = n.px / e_p;
    const double vy = n.py / e_p;
    const double vz = pz_p / e_p;
    const double v2 = vx * vx + vy * vy + vz * vz;
    if (v2 > d_V_max) d_V_max = v2;
  }
  for (const auto& n : b) {
    const double e = std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
    const double pz_p = g * (n.pz - beta * e);
    const double e_p  = g * (e - beta * n.pz);
    const double vx = n.px / e_p;
    const double vy = n.py / e_p;
    const double vz = pz_p / e_p;
    const double v2 = vx * vx + vy * vy + vz * vz;
    if (v2 > d_V_max) d_V_max = v2;
  }
  d_V_max = std::sqrt(d_V_max);

  const double dz_eff = g * dz;
  const double d_r = std::sqrt(dx * dx + dy * dy + dz_eff * dz_eff);

  return std::max(d_r, d_V_max * dmin / V0pass);
}

double GMSTClustering::CrossCMMetric(const std::vector<Node>& a, const std::vector<Node>& b,
                                     double dmin, double V0pass) const {
  // Центроиды групп
  double cxA = 0, cyA = 0, czA = 0;
  for (const auto& n : a) { cxA += n.x; cyA += n.y; czA += n.z; }
  const double invA = 1.0 / static_cast<double>(a.size());
  cxA *= invA; cyA *= invA; czA *= invA;

  double cxB = 0, cyB = 0, czB = 0;
  for (const auto& n : b) { cxB += n.x; cyB += n.y; czB += n.z; }
  const double invB = 1.0 / static_cast<double>(b.size());
  cxB *= invB; cyB *= invB; czB *= invB;

  const double dx = cxA - cxB;
  const double dy = cyA - cyB;
  const double dz = czA - czB;

  // Режим 0 — чисто пространственный
  if (metricMode == 0) {
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  // Если порог скорости выключен — ведём себя как режим 0
  if (V0pass <= 0.0) {
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  // 4-импульсы кластеров
  G4LorentzVector Pa4(0, 0, 0, 0), Pb4(0, 0, 0, 0);
  for (const auto& n : a) {
    const double e = std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
    Pa4 += G4LorentzVector(n.px, n.py, n.pz, e);
  }
  for (const auto& n : b) {
    const double e = std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
    Pb4 += G4LorentzVector(n.px, n.py, n.pz, e);
  }

  // Cross-CM per-particle cut: каждый нуклон проверяется в СЦМ другого кластера
  // (как в DCM BCOAV2 — кандидат проверяется в СО уже существующего кластера)
  double d_V_max = 0;
  const G4ThreeVector betaB = Pb4.boostVector();
  for (const auto& n : a) {
    const double e = std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
    G4LorentzVector Ln(n.px, n.py, n.pz, e);
    Ln.boost(-betaB);
    const double v2 = Ln.vect().mag2() / (Ln.e() * Ln.e());
    if (v2 > d_V_max) d_V_max = v2;
  }
  const G4ThreeVector betaA = Pa4.boostVector();
  for (const auto& n : b) {
    const double e = std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + n.mass * n.mass);
    G4LorentzVector Ln(n.px, n.py, n.pz, e);
    Ln.boost(-betaA);
    const double v2 = Ln.vect().mag2() / (Ln.e() * Ln.e());
    if (v2 > d_V_max) d_V_max = v2;
  }
  d_V_max = std::sqrt(d_V_max);

  // Пространственная часть — как в PairMetric (сжатие вдоль z в СЦМ объединения)
  double Ptz = Pa4.pz() + Pb4.pz();
  double Etot = Pa4.e() + Pb4.e();
  double beta = Ptz / Etot;
  if (beta >=  1.0) beta =  1.0 - 1e-12;
  if (beta <= -1.0) beta = -1.0 + 1e-12;
  const double g = 1.0 / std::sqrt(1.0 - beta * beta);

  const double dz_eff = g * dz;
  const double d_r = std::sqrt(dx * dx + dy * dy + dz_eff * dz_eff);

  return std::max(d_r, d_V_max * dmin / V0pass);
}

std::vector<std::vector<G4int> > GMSTClustering::MergeClusters(
    const std::vector<std::vector<G4int> >& clusters,
    const aamcc::NucleonVector* nucleons,
    double d0, double V0pass) const {
  const G4int K = static_cast<G4int>(clusters.size());
  if (K <= 1) return clusters;

  // Узлы-кластеры: вектор нуклонов каждого кластера (для per-constituent метрики)
  std::vector<std::vector<Node> > clusterNodes(K);
  for (G4int c = 0; c < K; ++c) {
    clusterNodes[c].reserve(clusters[c].size());
    for (G4int idx : clusters[c]) {
      const aamcc::Nucleon* n = &nucleons->at(idx);
      clusterNodes[c].push_back({ n->GetX(), n->GetY(), n->GetZ(),
                                  n->GetPx(), n->GetPy(), n->GetPz(),
                                  nucleonAverMass, 1 });
    }
  }

  // Рёбра с индивидуальным порогом dmin_eff = d0 * (A_i + A_j)^(1/3)
  // DCM-подобная метрика: per-constituent |V_k'| < V0 в СЦМ объединения
  std::vector<std::pair<double, iPair> > edges;
  std::vector<double> edgeDmin;
  edges.reserve(static_cast<size_t>(K) * (K - 1) / 2);
  edgeDmin.reserve(static_cast<size_t>(K) * (K - 1) / 2);
  for (G4int i = 0; i < K; ++i) {
    for (G4int j = i + 1; j < K; ++j) {
      const double dm = d0 * std::pow(static_cast<double>(clusters[i].size() + clusters[j].size()), 1.0 / 3.0);
      edges.push_back({PairMetric(clusterNodes[i], clusterNodes[j], dm, V0pass), {i, j}});
      edgeDmin.push_back(dm);
    }
  }
  std::sort(edges.begin(), edges.end());

  // MST + отсечение по индивидуальному порогу (та же логика, что в Graph)
  DisjointSets ds(K);
  std::vector<std::vector<G4int> > adj(K);
  for (size_t e = 0; e < edges.size(); ++e) {
    const G4int u = edges[e].second.first;
    const G4int v = edges[e].second.second;
    const G4int su = ds.find(u);
    const G4int sv = ds.find(v);
    if (su != sv) {
      if (edges[e].first < edgeDmin[e]) { adj[u].push_back(v); adj[v].push_back(u); }
      ds.merge(su, sv);
    }
  }

  // Компоненты связности (DFS) -> объединённые кластеры
  std::vector<std::vector<G4int> > merged;
  std::vector<G4bool> visited(K, false);
  for (G4int v = 0; v < K; ++v) {
    if (visited[v]) continue;
    std::vector<G4int> comp;
    std::stack<G4int> st;
    st.push(v); visited[v] = true;
    while (!st.empty()) {
      const G4int x = st.top(); st.pop();
      comp.push_back(x);
      for (G4int y : adj[x]) { if (!visited[y]) { visited[y] = true; st.push(y); } }
    }
    std::vector<G4int> nucs;
    for (G4int c : comp) for (G4int idx : clusters[c]) nucs.push_back(idx);
    merged.push_back(nucs);
  }
  return merged;
}

void GMSTClustering::SetCDExEn(G4double Ex, G4int A) {
  // ==== Установка CritDist для режимов 1 и 2 ====
  // В режимах pair-CM CritDist фиксируется на d0 (пороги задаются dmin и V0)
  if (metricMode != 0) {
    CritDist = d0;
    return;
  }

  if (Ex / G4double(A) < 2.17 * MeV) {
    CritDist = d0;
  } else {
    G4double ex = Ex / G4double(A);
    G4double dep = std::exp(-std::pow((ex / b_opt), a_opt)) * c_opt + d_opt;
    CritDist = d0 * std::pow(dep, 1. / 3.);
  }
  if (alphaPow > 10) {
    CritDist = d0;
  }
}

std::vector<G4FragmentVector> GMSTClustering::GetClusters(aamcc::NucleonVector *nucleons_in, G4double ExA, G4double ExB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB) {
  G4int A = 0;
  G4int Z = 0;
  G4int Ab = 0;
  G4int Zb = 0;
  nucleonVector = nucleons_in;
  SpecAa = nucleonVector->GetA("A");
  SpecAb = nucleonVector->GetA("B");

  auto *nucleons = new aamcc::NucleonVector();    // nucleons from Side A
  auto *nucleons_B = new aamcc::NucleonVector();  // nucleons from Side B

  for (G4int iArray = 0; iArray < nucleons_in->size(); iArray++) {
    aamcc::Nucleon *nucleon = &(nucleons_in->at(iArray));
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "A") {
      A += 1;
      nucleons->push_back(nucleons_in->at(iArray));
    }
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "A" && nucleon->isospin == 1) {
      Z += 1;
    }
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "B") {
      Ab += 1;
      nucleons_B->push_back(nucleons_in->at(iArray));
    }
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "B" && nucleon->isospin == 1) {
      Zb += 1;
    }
  }
  // Graph для side A
  Graph g;

  // Graph для side B
  Graph g_B;

  vector<vector<G4int> > clusters;
  vector<vector<G4int> > clusters_B;

  const int npass = (metricMode == 2) ? 2 : 1;

  if (npass == 1) {
    if (ExA > 0 && ExB > 0) this->SetCDExEn(ExA, A);
    CritDistA = this->CritDist;
    g = this->ClusterToGraph(nucleons, V0);
    clusters = g.AdvancedKruskalMST(this->CritDist);

    if (ExA > 0 && ExB > 0) this->SetCDExEn(ExB, Ab);
    g_B = this->ClusterToGraph(nucleons_B, V0);
    clusters_B = g_B.AdvancedKruskalMST(this->CritDist);
  } else {
    // mode 2: 2 иерархических прохода, постоянный V0
    // проход 1: V0, dmin = d0 (CritDist фиксирован)
    if (ExA > 0 && ExB > 0) this->SetCDExEn(ExA, A);
    CritDistA = this->CritDist;
    {
      Graph g1 = this->ClusterToGraph(nucleons, V0);
      clusters = g1.AdvancedKruskalMST(this->CritDist);
    }
    // проход 2: V0, dmin_eff = d0*(A_i+A_j)^(1/3)
    clusters = MergeClusters(clusters, nucleons, d0, V0);

    // side B
    if (ExA > 0 && ExB > 0) this->SetCDExEn(ExB, Ab);
    {
      Graph g1b = this->ClusterToGraph(nucleons_B, V0);
      clusters_B = g1b.AdvancedKruskalMST(this->CritDist);
    }
    clusters_B = MergeClusters(clusters_B, nucleons_B, d0, V0);
  }

  std::vector<G4FragmentVector> outClusters;
  G4FragmentVector output_vector_A;
  G4FragmentVector output_vector_B;

  std::vector<int> rmapsA;  // repulsion maps with fragments for size A
  std::vector<int> rmapsB;  // repulsion maps with fragments for size B
  int rcountA = 0;
  int rcountB = 0;
  aamcc::NucleonVector rnucsA;  // protons for size A
  aamcc::NucleonVector rnucsB;  // protons for size B

  // Filling a Cluster Vector (side A)
  for (G4int i = 0; i < clusters.size(); ++i) {
    G4int Z_clust = 0;
    G4int A_clust = 0;
    for (G4int j = 0; j < clusters[i].size(); ++j) {
      aamcc::Nucleon *nucleon = &(nucleons->at((clusters[i])[j]));

      if (nucleon->isospin == 1) {
        rnucsA.push_back(*nucleon);
        rmapsA.push_back(rcountA);

        Z_clust += 1;
      }
      A_clust += 1;
    }

    G4Fragment *frag = new G4Fragment();
    frag->SetA(A_clust);
    frag->SetZ(Z_clust);
    output_vector_A.push_back(frag);

    ++rcountA;
  }

  // Filling a Cluster Vector (side B)
  for (G4int i = 0; i < clusters_B.size(); ++i) {
    G4int Z_clust = 0;
    G4int A_clust = 0;
    for (G4int j = 0; j < clusters_B[i].size(); ++j) {
      aamcc::Nucleon *nucleon = &(nucleons_B->at((clusters_B[i])[j]));

      if (nucleon->isospin == 1) {
        rnucsB.push_back(*nucleon);
        rmapsB.push_back(rcountB);

        Z_clust += 1;
      }
      A_clust += 1;
    }

    G4Fragment *frag = new G4Fragment();
    frag->SetA(A_clust);
    frag->SetZ(Z_clust);
    output_vector_B.push_back(frag);

    ++rcountB;
  }
  outClusters.push_back(output_vector_A);
  outClusters.push_back(output_vector_B);

  delete nucleons;
  delete nucleons_B;

  return CalculateMomentum(outClusters, ExA, ExB, boostA, boostB, rnucsA, rnucsB, rmapsA, rmapsB);
};

std::vector<G4FragmentVector> GMSTClustering::CalculateMomentum(std::vector<G4FragmentVector> noMomClusters, G4double ExEnA, G4double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, aamcc::NucleonVector rnucsA, aamcc::NucleonVector rnucsB, std::vector<int> rmapsA, std::vector<int> rmapsB) {
  std::vector<G4FragmentVector> momClusters = noMomClusters;
  std::vector<G4double> MstMassVector_A;
  MstMassVector_A.reserve(noMomClusters.at(0).size());

  std::vector<G4double> MstMassVector_B;
  MstMassVector_B.reserve(noMomClusters.at(1).size());

  G4double SumMassMst = 0;
  G4double SumMassMstEx = 0;

  if (!noMomClusters.at(0).empty()) {
    for (G4int i = 0; i < noMomClusters.at(0).size(); ++i) {
      G4int clfrag_A = (noMomClusters.at(0)[i])->GetA();
      G4int clfrag_Z = (noMomClusters.at(0)[i])->GetZ();
      // Excitation energy computing
      G4double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnA * G4double(clfrag_A) / G4double(SpecAa);
      }

      G4double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

      SumMassMstEx += NuclearMass;

      MstMassVector_A.push_back(NuclearMass);
    }

    G4double PrefragmentMass_A = SumMassMst + ExEnA;

    std::vector<G4LorentzVector *>* momentumVectorA;
    if (PrefragmentMass_A < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_A += 1e-5 * MeV;
    }
    momentumVectorA = phaseSpaceDecay.Decay(PrefragmentMass_A, MstMassVector_A);

    for (int I = 0; I < momClusters.at(0).size(); ++I) {
      momClusters.at(0).at(I)->SetMomentum(*momentumVectorA->at(I));
    }

    momentumVectorA->clear();

    momClusters.at(0) = RepulsionStage::CalculateRepulsion(momClusters.at(0), rnucsA, rmapsA);
    for (int I = 0; I < momClusters.at(0).size(); ++I) {
      G4LorentzVector momentum = momClusters.at(0).at(I)->GetMomentum();
      momentum.boost(boostA);
      momClusters.at(0).at(I)->SetMomentum(momentum);
    }
  }
  // side B
  SumMassMst = 0;
  SumMassMstEx = 0;
  if (!noMomClusters.at(1).empty()) {
    for (G4int i = 0; i < noMomClusters.at(1).size(); ++i) {
      G4int clfrag_A = (noMomClusters.at(1)[i])->GetA();
      G4int clfrag_Z = (noMomClusters.at(1)[i])->GetZ();
      // Excitation energy computing
      G4double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnB * G4double(clfrag_A) / G4double(SpecAb);
      }

      G4double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

      SumMassMstEx += NuclearMass;

      MstMassVector_B.push_back(NuclearMass);
    }

    G4double PrefragmentMass_B = SumMassMst + ExEnB;
    std::vector<G4LorentzVector *> *momentumVectorB;
    if (PrefragmentMass_B < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_B += 1e-5 * MeV;
    }
    momentumVectorB = phaseSpaceDecay.Decay(PrefragmentMass_B, MstMassVector_B);

    for (int I = 0; I < momClusters.at(1).size(); ++I) {
      momClusters.at(1).at(I)->SetMomentum(*momentumVectorB->at(I));
    }

    momentumVectorB->clear();

    momClusters.at(1) = RepulsionStage::CalculateRepulsion(momClusters.at(1), rnucsB, rmapsB);
    for (int I = 0; I < momClusters.at(1).size(); ++I) {
      G4LorentzVector momentum = momClusters.at(1).at(I)->GetMomentum();
      momentum.boost(boostB);
      momClusters.at(1).at(I)->SetMomentum(momentum);
    }
  }

  MstMassVector_A.clear();
  MstMassVector_B.clear();
  noMomClusters.clear();

  return momClusters;
};

Graph::Graph(G4int V_in, G4int E_in) {
  this->V = V_in;
  this->E = E_in;
  adj.assign((size_t)V_in, list<G4int>());
}

Graph::Graph() {
  this->V = 0;
  this->E = 0;
}

Graph::Graph(const Graph& other) {
  V = other.V;
  E = other.E;
  edges = other.edges;
  adj = other.adj;
}

Graph& Graph::operator=(const Graph& other) {
  if (this != &other) {
    V = other.V;
    E = other.E;
    edges = other.edges;
    adj = other.adj;
  }
  return *this;
}

Graph::Graph(Graph&& other) noexcept {
  V = other.V; E = other.E;
  edges = std::move(other.edges);
  adj = std::move(other.adj);
  other.V = 0; other.E = 0;
}

Graph& Graph::operator=(Graph&& other) noexcept {
  if (this != &other) {
    V = other.V; E = other.E;
    edges = std::move(other.edges);
    adj = std::move(other.adj);
    other.V = 0; other.E = 0;
  }
  return *this;
}

Graph::~Graph() {}

void Graph::addEdge(G4int u, G4int v, G4double w) { edges.push_back({w, {u, v}}); }

vector<vector<G4int> > Graph::AdvancedKruskalMST(G4double CD_in) {
  // Sort edges in increasing order on basis of cost
  sort(edges.begin(), edges.end());

  // Create disjoint sets
  DisjointSets ds(V);

  // Iterate through all sorted edges
  vector<pair<G4double, iPair> >::iterator it;
  for (it = edges.begin(); it != edges.end(); it++) {
    G4int u = it->second.first;
    G4int v = it->second.second;

    G4int set_u = ds.find(u);
    G4int set_v = ds.find(v);

    // Check if the selected edge is creating
    // a cycle or not (Cycle is created if u
    // and v belong to same set)
    if (set_u != set_v) {
      // Current edge will be in the MST
      if (it->first < CD_in) {
        this->addConn(u, v);
      }
      // Merge two sets
      ds.merge(set_u, set_v);
    }
  }

  return this->connectedComponents();
}

void Graph::addConn(G4int v, G4int w) {
  if (V <= 0) return;
  if (v < 0 || w < 0) return;
  if (v >= V || w >= V) return;
  adj[(size_t)v].push_back(w);
  adj[(size_t)w].push_back(v);
}

vector<vector<G4int> > Graph::connectedComponents() {
  vector<vector<G4int> > clusters;
  vector<vector<G4int> > clusters_B;
  // Mark all the vertices as not visited
  vector<char> visited((size_t)V, 0);
  for (G4int v = 0; v < V; v++) {
    if (!visited[(size_t)v]) {
      // vector to fill with nucleon's numbers from one cluster
      vector<G4int> clust_inside;
      // print all reachable vertices from v
      DFSUtil(v, (G4bool*)visited.data(), &clust_inside);
      clusters.push_back(clust_inside);
    }
  }
  return clusters;
}

void Graph::DFSUtil(G4int v, G4bool visited[], vector<G4int> *clust_inside) {
  // Mark the current node as visited and print it
  visited[(size_t)v] = true;
  (*clust_inside).push_back(v);

  // Recur for all the vertices adjacent to this vertex
  list<G4int>::iterator i;
  for (i = adj[(size_t)v].begin(); i != adj[(size_t)v].end(); ++i) {
    if (*i >= 0 && *i < V && !visited[(size_t)(*i)]) DFSUtil(*i, visited, &(*clust_inside));
  }
}

DisjointSets::DisjointSets(G4int nsz) {
  // Allocate memory
  this->n = nsz;
  parent = new G4int[n + 1];
  rnk = new G4int[n + 1];

  // Initially, all vertices are in different sets and have rank 0.
  for (G4int i = 0; i <= n; i++) {
    rnk[i] = 0;
    // every element is parent of itself
    parent[i] = i;
  }
}

DisjointSets::~DisjointSets() {
  delete[] parent;
  delete[] rnk;
}

G4int DisjointSets::find(G4int u) {
  /* Make the parent of the nodes in the path
     from u--> parent[u] point to parent[u] */
  if (u != parent[u]) parent[u] = find(parent[u]);
  return parent[u];
}

void DisjointSets::merge(G4int x, G4int y) {
  x = find(x), y = find(y);

  /* Make tree with smaller height
     a subtree of the other tree  */
  if (rnk[x] > rnk[y])
    parent[y] = x;
  else  // If rnk[x] <= rnk[y]
    parent[x] = y;

  if (rnk[x] == rnk[y]) rnk[y]++;
}