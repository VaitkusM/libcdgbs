#include "libcdgbs/SurfGBS.hpp"
#include <OpenMesh/Core/IO/MeshIO.hh>            // write_mesh()
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>  // registers only the OBJ writer

// force static registration of the OBJ writer
static OpenMesh::IO::_OBJWriter_ _objWriterRegistration_;

using namespace libcdgbs;

bool SurfGBS::writeOBJ(const Mesh& mesh, const std::string& filename)
{
  return OpenMesh::IO::write_mesh(mesh, filename);
}