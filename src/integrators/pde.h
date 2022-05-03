#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef INTEGRATORS_PDE_H
#define INTEGRATORS_PDE_H

#include "pbrt.h"
#include "integrator.h"
#include "scene.h"
#include "extension/pde/pdes.h"
#include "extension/config.h"
using namespace pbrt;

template <size_t N>
class Solver;

namespace extension {
    class PDEIntegrator : public SamplerIntegrator {
        public:
            // PDEIntegrator Public Methods
            PDEIntegrator(std::shared_ptr<const Camera> camera,
                std::shared_ptr<Sampler> sampler,
                const Bounds2i &pixelBounds, std::shared_ptr<Primitive> &boundary)
                : SamplerIntegrator(camera, sampler, pixelBounds), boundary(boundary),
                config{"C:/Personal/Project/experiment/pbrt-v3/scenes/config.txt"} {
                maxDist = 1.0 / (2 * config.numCol);
            }
            Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;
            void Preprocess(const Scene &scene, Sampler &sampler);
            void Postprocess(const Scene &scene, Sampler &sampler);
        private:
            std::shared_ptr<Primitive> boundary;
            Solver<3> *solver = nullptr;
            mutable SolverData<3> solverData;
            mutable PDEData<3> pdeData;
            mutable int numPointEvaluated = 0;
            Config config;
            double maxDist;
    };

	PDEIntegrator *CreatePDEIntegrator(
		const ParamSet &params, std::shared_ptr<Sampler> sampler,
		std::shared_ptr<const Camera> camera, std::shared_ptr<Primitive> boundary);

} // namespace extension
#endif  // INTEGRATORS_PDE_H