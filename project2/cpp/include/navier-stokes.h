
extern VectorField u;

extern FaceVectorField uf;
extern FaceVectorField mu;

extern ScalarField p;

static void advection_term (int istep, double t, double dt);


static void viscous_term (int istep, double t, double dt);


static void pressure_correction (int istep, double t, double dt);

