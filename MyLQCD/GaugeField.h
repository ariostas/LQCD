template <class T>
class GaugeField { 
   std::vector<T> links;
   size_t Nx, Ny, Nz, Nt;
public:
    GaugeField(size_t Nx, size_t Ny, size_t Nz, size_t Nt) 
        : links(Nx*Ny*Nz*Nt*4),
        Nx(Nx), Ny(Ny), Nz(Nz), Nt(Nt)
    {}

    T &operator()(size_t nx, size_t ny, size_t nz, size_t nt, size_t dir) {
        return links.at( (nx%Nx)*(Ny*Nz*Nt*4) + (ny%Ny)*(Nz*Nt*4) + (nz%Nz)*(Nt*4) + (nt%Nt)*(4) + (dir%4) );
    }
};