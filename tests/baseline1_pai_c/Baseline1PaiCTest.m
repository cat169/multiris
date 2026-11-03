classdef Baseline1PaiCTest < matlab.unittest.TestCase
    properties (TestParameter)
        pattern={"111","110","101","011"}
    end
    methods (TestClassSetup)
        function paths(test)
            repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
            test.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo,'src/baselines/multiris3d')));
        end
    end
    methods (Test)
        function residualScaling(test,pattern)
            [cfg,ch,x]=PaiCFixture.data(pattern);
            stream=RandStream('mt19937ar','Seed',71303);x=x+.1*randn(stream,6,1);
            original=PaiCFixture.secondsResidual(x,ch,cfg);
            [scaled,~]=baseline1_distance_scaled_residual(x,ch.delayEstimate,ch.xiEstimate,ch.zetaEstimate,ch.activeRisIndices,cfg);
            test.verifyEqual(scaled([1:3:end 2:3:end]),original([1:3:end 2:3:end]),'AbsTol',0);
            test.verifyEqual(scaled(3:3:end),cfg.constants.c*original(3:3:end),'AbsTol',1e-13);
            test.verifyEqual(sum(scaled.^2),PaiCFixture.independentCost(x,ch,cfg),'AbsTol',1e-12);
        end
        function analyticJacobian(test,pattern)
            [absolute,relative]=PaiCFixture.finiteDifference(pattern,false);
            test.verifyLessThan(absolute,1e-7);test.verifyLessThan(relative,1e-7);
        end
        function refinerCosts(test,pattern)
            [cfg,ch,x]=PaiCFixture.data(pattern);
            result=refine_position_mle_eq57_distance_scaled(ch.delayEstimate,ch.xiEstimate,ch.zetaEstimate, ...
                ch.activeRisIndices,x(1:3),x(4:6),cfg);
            final=[result.staticPositionEstimate;result.mobilePositionEstimate];
            test.verifyEqual(result.initialCost,PaiCFixture.independentCost(x,ch,cfg),'AbsTol',1e-11);
            test.verifyEqual(result.finalCost,PaiCFixture.independentCost(final,ch,cfg),'AbsTol',1e-11);
            test.verifyEqual(result.residualOrder,["xi";"zeta";"distance_metres"]);
        end
        function normalPathObjective(test)
            result=PaiCFixture.entry(false);
            [a,b,~,selected]=PaiCFixture.candidateChecks(result);
            test.verifyFalse(result.fallbackTriggered);test.verifyTrue(result.finite);
            test.verifyEqual(result.candidateCount,1);test.verifyTrue(selected);
            test.verifyLessThan(max(a),1e-11);test.verifyLessThan(max(b),1e-11);
            test.verifyEqual(result.candidates.spatialObjective+result.candidates.delayObjective,result.finalObjective,'AbsTol',1e-12);
        end
        function fallbackObjectiveAndSelection(test)
            result=PaiCFixture.entry(true);
            [a,b,sorted,selected]=PaiCFixture.candidateChecks(result);
            test.verifyTrue(result.fallbackTriggered);test.verifyTrue(result.finite);
            test.verifyGreaterThan(result.candidateCount,1);test.verifyLessThanOrEqual(result.candidateCount,8);
            test.verifyTrue(sorted);test.verifyTrue(selected);
            test.verifyLessThan(max(a),1e-11);test.verifyLessThan(max(b),1e-11);
        end
        function nearCenterReportedSeparately(test)
            [absolute,relative]=PaiCFixture.finiteDifference("111",true);
            test.verifyTrue(isfinite(absolute)&&isfinite(relative));
        end
        function exactCenterRejected(test)
            [cfg,ch,x]=PaiCFixture.data("111");x(1:3)=cfg.geometry.ris_centers(:,1);
            test.verifyError(@()baseline1_distance_scaled_residual(x,ch.delayEstimate,ch.xiEstimate,ch.zetaEstimate,ch.activeRisIndices,cfg), ...
                'refine_position_mle_eq57_distance_scaled:PositionAtRisCenter');
        end
        function originalAlgorithmsUnchanged(test)
            test.verifyTrue(PaiCFixture.protectedHashes());
        end
    end
end
