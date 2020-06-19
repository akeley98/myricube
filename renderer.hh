#ifndef MYRICUBE_RENDERER_HH_
#define MYRICUBE_RENDERER_HH_

namespace myricube {

class TextureStore;
class MeshStore;

TextureStore* new_texture_store();
MeshStore* new_mesh_store();
void delete_texture_store(TextureStore*);
void delete_mesh_store(MeshStore*);

} // end namespace
#endif /* !MYRICUBE_RENDERER_HH_ */
