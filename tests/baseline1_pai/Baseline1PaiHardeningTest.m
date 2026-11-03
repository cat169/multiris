classdef Baseline1PaiHardeningTest < matlab.unittest.TestCase
    properties
        Fixture
    end
    methods (TestMethodSetup)
        function isolateSources(testCase)
            testCase.Fixture=PaiHardeningFixture();
            testCase.addTeardown(@() testCase.Fixture.close());
        end
    end
    methods (Test)
        function rejectNonstructChannel(t)
            f=t.Fixture;
            t.verifyError(@() infer_baseline1_frame_projection_assisted([],1,[],f.Configuration,42),'baseline1PAI:InvalidExternalChannel');
        end
        function rejectMissingField(t)
            f=t.Fixture; f.Channel=rmfield(f.Channel,'xiEstimate');
            t.verifyError(@() f.call(),'baseline1PAI:InvalidExternalChannel');
        end
        function rejectSupportOrder(t)
            f=t.Fixture; f.Channel.activeRisIndices=[2;1];
            t.verifyError(@() f.call(),'baseline1PAI:ExternalSupportMismatch');
        end
        function rejectWrongVectorLength(t)
            f=t.Fixture; f.Channel.xiEstimate=1;
            t.verifyError(@() f.call(),'baseline1PAI:ExternalChannelSize');
        end
        function rejectMatrixEstimate(t)
            f=t.Fixture; f.Channel.xiEstimate=ones(2);
            t.verifyError(@() f.call(),'baseline1PAI:ExternalChannelSize');
        end
        function rejectNaN(t)
            f=t.Fixture; f.Channel.delayEstimate(1)=nan;
            t.verifyError(@() f.call(),'baseline1PAI:ExternalChannelNonfinite');
        end
        function rejectInf(t)
            f=t.Fixture; f.Channel.zetaEstimate(1)=inf;
            t.verifyError(@() f.call(),'baseline1PAI:ExternalChannelNonfinite');
        end
        function rejectComplex(t)
            f=t.Fixture; f.Channel.xiEstimate(1)=1i;
            t.verifyError(@() f.call(),'baseline1PAI:ExternalChannelNonfinite');
        end
        function rejectNoiseMismatch(t)
            f=t.Fixture; f.Channel.separatedNoiseVariance=2;
            t.verifyError(@() f.call(),'baseline1PAI:ExternalNoiseMismatch');
        end
        function skipExactRisCenter(t)
            f=t.Fixture; f.Configuration.paperBaselineA.lowerBound=[0;0;0];
            r=f.call();
            t.verifyTrue(r.finite); t.verifyTrue(r.initializationDiagnostics.txAtRisCenter(1));
            t.verifyEqual(nnz(r.initializationDiagnostics.txAtRisCenter),1);
        end
        function continueAfterNaNPosition(t)
            f=t.Fixture; f.setMode('singleNonfinite'); r=f.call();
            t.verifyTrue(r.finite); t.verifyEqual(r.selectedIndex,2);
            t.verifyFalse(r.candidates(1).valid); t.verifyTrue(all(isnan(r.candidates(1).estimate(1:3))));
            t.verifyEqual(r.candidates(1).exceptionIdentifier,'baseline1PAI:NonfiniteCandidate');
        end
        function continueAfterInfFinalObjective(t)
            f=t.Fixture; f.setMode('singleInfObjective'); r=f.call();
            t.verifyTrue(r.finite); t.verifyFalse(r.candidates(1).valid);
            t.verifyEqual(r.candidates(1).finalObjective,r.candidates(1).formal.initialCost,'AbsTol',1e-12);
            t.verifyEqual(r.candidates(1).formal.finalCost,inf);
        end
        function continueAfterInfInitialObjective(t)
            f=t.Fixture; f.setMode('nonfiniteInitial'); r=f.call();
            t.verifyTrue(r.finite); t.verifyFalse(r.candidates(1).valid);
            t.verifyEqual(r.candidates(1).formal.initialCost,inf);
        end
        function allCandidatesInvalidReturnsFailure(t)
            f=t.Fixture; f.setMode('allNonfinite'); r=f.call();
            t.verifyFalse(r.finite); t.verifyTrue(r.fallbackTriggered); t.verifyEqual(r.selectedIndex,0);
            t.verifyEqual(r.exceptionIdentifier,'baseline1PAI:NoValidRefinement');
            t.verifyTrue(all(isnan([r.staticPositionEstimate;r.mobilePositionEstimate])));
            t.verifyEmpty(fieldnames(r.coarse)); t.verifyEmpty(fieldnames(r.refined));
        end
        function allGridPointsInvalidReturnsFailure(t)
            f=t.Fixture; f.Configuration.paperBaselineA.lowerBound=[0;0;0]; f.Configuration.paperBaselineA.upperBound=[0;0;0]; r=f.call();
            t.verifyFalse(r.finite); t.verifyEqual(r.candidateCount,0);
            t.verifyEqual(r.exceptionIdentifier,'baseline1PAI:NoValidRefinement');
        end
        function initialMismatchIsFatal(t)
            f=t.Fixture; f.setMode('initialMismatch');
            t.verifyError(@() f.call(),'fallback:InitialObjectiveMismatch');
        end
        function finalMismatchIsFatal(t)
            f=t.Fixture; f.setMode('finalMismatch');
            t.verifyError(@() f.call(),'fallback:FinalObjectiveMismatch');
        end
        function dimensionErrorIsFatal(t)
            f=t.Fixture; f.setMode('dimensionError');
            t.verifyError(@() f.call(),'MATLAB:dimagree');
        end
        function unknownErrorIsFatal(t)
            f=t.Fixture; f.setMode('unknownError');
            t.verifyError(@() f.call(),'paiTest:UnknownProgramError');
        end
        function knownNumericalExceptionContinues(t)
            f=t.Fixture; f.setMode('numericException'); r=f.call();
            t.verifyTrue(r.finite); t.verifyFalse(r.candidates(1).valid);
            t.verifyEqual(r.candidates(1).exceptionIdentifier,'refine_position_mle_eq57_identity:PositionAtRisCenter');
        end
        function nonconvergenceDoesNotFilterFiniteCandidate(t)
            f=t.Fixture; f.setMode('finiteNotConverged'); r=f.call();
            t.verifyTrue(r.finite); t.verifyTrue(r.qualityFailed); t.verifyEqual(r.refined.exitflag,0);
        end
        function fourArgumentPathEstimatesMleOnce(t)
            f=t.Fixture; f.setMode('normal'); r=f.callDefault(); s=f.state();
            t.verifyFalse(r.fallbackTriggered); t.verifyEqual(s.channelCalls,1); t.verifyEqual(s.coarseCalls,1); t.verifyEqual(s.refineCalls,1);
            t.verifyEqual(r.extraPositionSolves,0);
        end
        function threeReturnCasesHaveSameTopLevelFields(t)
            f=t.Fixture; f.setMode('normal'); normal=f.call(); f.setMode('finite'); fallback=f.call(); f.setMode('allNonfinite'); failed=f.call();
            t.verifyEqual(fieldnames(normal),fieldnames(fallback)); t.verifyEqual(fieldnames(fallback),fieldnames(failed));
        end
    end
end
