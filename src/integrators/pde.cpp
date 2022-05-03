#include "integrators/pde.h"
#include "extension/boundary/pbrtboundary.h"
#include "extension/methods/wos.h"

#include "interaction.h"
#include "camera.h"
#include "film.h"
#include "paramset.h"

using namespace pbrt;
using namespace MyTypes;

namespace extension {
    // PDEIntegrator Method Definitions

    void PDEIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
        TaylorVortexPDE3 *pde = new TaylorVortexPDE3();
        PBRTBoundary *bb = new PBRTBoundary([](arrayd<3> const &) {return 0;}, boundary, config.maxError);
        std::random_device rd;
        switch (config.solverType) {
        case SolverTypes::WOS:
            solver = new WOS3d(pde, bb, std::mt19937(rd()));
            break;
        default:
            throw exception("Unknown solver type");
        }
        if (config.refName != "normal") {
            config.getSol(pdeData.streamRef);
        }
        solverData.streamError.resize(config.stopAt);
    }

    void PDEIntegrator::Postprocess(const Scene &scene, Sampler &sampler) {
        config.saveSol(solverData.streamAcc);
    }

    Spectrum PDEIntegrator::Li(const RayDifferential &ray, const Scene &scene,
        Sampler &sampler, MemoryArena &arena, int depth) const {
        Spectrum L(0.);
        // Find closest ray intersection or return background radiance
        SurfaceInteraction isect;
        if (!scene.Intersect(ray, &isect)) {
            for (const auto &light : scene.lights) L += light->Le(ray);
            return L;
        }
        const Normal3f &n = isect.shading.n;
        Vector3f wo = isect.wo;
        // Compute scattering functions for surface interaction
        isect.ComputeScatteringFunctions(ray, arena);
        if (!isect.bsdf)
            return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
        Point3f x = isect.p;
        int localCoords[3]{ x.x * config.numCol, x.y * config.numRow, x.z * config.numZ };
        if (localCoords[0] > config.numCol - 1) localCoords[0] = config.numCol - 1;
        if (localCoords[1] > config.numRow - 1) localCoords[1] = config.numRow - 1;
        if (localCoords[2] > config.numZ - 1) localCoords[2] = config.numZ - 1;
        int i = localCoords[0] + localCoords[1] * config.numRow + localCoords[2] * config.numCol * config.numRow;
        CHECK_LT(i, config.numCells);
        double sol = 0;
        pdeData.streamLock.lock();
        if (solverData.streamAcc.count(i) == 0) {
            pdeData.computePoints[i] = { {x.x,x.y,x.z} };
            if (config.refName == "normal") {
                pdeData.streamRef[i] = solver->truth(pdeData.computePoints[i], 0);
            }
            if (boundary->minDistanceFromPoint(x) > config.maxError ) {
                for (int s =  1; s <= config.stopAt; s++) {
                    solverData.streamAcc[i] += solver->eval(pdeData.computePoints[i], 0);
                    solverData.streamError[s - 1] += pdeData.streamRef[i] - solverData.streamAcc[i] / s;
                }
                solverData.streamAcc[i] /= config.stopAt;
            }
            else {
                solverData.streamAcc[i] = solver->cond(pdeData.computePoints[i]);
            }
            numPointEvaluated++;
        }
        pdeData.streamLock.unlock();
        sol = solverData.streamAcc[i];
        Float color[3] = { 0, 0, 0 };
        if (sol > 0) color[0] = sol;
        else color[2] = -sol;
        // Add contribution of each light source
        for (const auto &light : scene.lights) {
            Vector3f wi;
            Float pdf;
            VisibilityTester visibility;
            Spectrum Li =
                light->Sample_Li(isect, sampler.Get2D(), &wi, &pdf, &visibility);
            if (Li.IsBlack() || pdf == 0) continue;
            Spectrum val = Spectrum::FromRGB(color);
            if (visibility.Unoccluded(scene)){
                if (boundary->minDistanceFromPoint(x) <= config.maxError) {
                    ray.tMax = Infinity;
                    boundary->Intersect(ray, &isect);
                    isect.ComputeScatteringFunctions(ray, arena);
                    Spectrum f = isect.bsdf->f(wo, wi);
                    L += f * Li * AbsDot(wi, n) / pdf;
                }
                else L += val * InvPi * Li * AbsDot(wi, n) / pdf;
            }
        }
        return L;
    }

    PDEIntegrator *CreatePDEIntegrator(
        const ParamSet &params, std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera, std::shared_ptr<Primitive> boundary) {
        int np;
        const int *pb = params.FindInt("pixelbounds", &np);
        Bounds2i pixelBounds = camera->film->GetSampleBounds();
        if (pb) {
            if (np != 4)
                Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                    np);
            else {
                pixelBounds = Intersect(pixelBounds,
                    Bounds2i{ {pb[0], pb[2]}, {pb[1], pb[3]} });
                if (pixelBounds.Area() == 0)
                    Error("Degenerate \"pixelbounds\" specified.");
            }
        }
        return new PDEIntegrator(camera, sampler, pixelBounds, boundary);
    }

} // namespace extension