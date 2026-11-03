classdef Baseline2S4ContractTest < matlab.unittest.TestCase
    properties
        Fixture
    end
    methods(TestMethodSetup)
        function setup(testCase)
            testCase.Fixture=S4ContractFixture();
            testCase.addTeardown(@()testCase.Fixture.close());
        end
    end
    methods(Test)
        function allFourValid(testCase)
            r=testCase.Fixture.call();
            testCase.verifyTrue(r.finite);
            testCase.verifyEqual(r.candidateStatus,"allCandidatesValid");
            testCase.verifyEqual(r.selectedIndex,1);
            testCase.verifyEqual(r.callCounts.map,4);
        end
        function recoverablePartialFailure(testCase)
            testCase.Fixture.setMode('partial');r=testCase.Fixture.call();
            testCase.verifyTrue(r.finite);
            testCase.verifyEqual(r.candidateStatus,"partialCandidateFailure");
            testCase.verifyEqual(r.candidates(2).exceptionIdentifier,"refine_position_map:PositionAtRisCenter");
        end
        function allCandidatesFailed(testCase)
            testCase.Fixture.setMode('allfail');r=testCase.Fixture.call();
            testCase.verifyFalse(r.finite);
            testCase.verifyEqual(r.selectedIndex,0);
            testCase.verifyEqual(r.candidateStatus,"allCandidatesFailed");
            testCase.verifyEqual(r.exceptionIdentifier,"baseline2S4:NoValidCandidate");
            testCase.verifyTrue(all(isnan(r.positionEstimate)));
            testCase.verifyEmpty(fieldnames(r.map));
        end
        function originalSrFailureStillAllowsOtherStarts(testCase)
            testCase.Fixture.setMode('firstfail');r=testCase.Fixture.call();
            testCase.verifyTrue(r.finite);
            testCase.verifyEqual(r.selectedIndex,2);
            testCase.verifyEqual(r.candidateStatus,"partialCandidateFailure");
            testCase.verifyEqual(r.srOriginalException.identifier,"refine_position_map:PositionAtRisCenter");
        end
        function failedSwapPreservesSrAndOtherStarts(testCase)
            testCase.Fixture.setMode('swapfail');r=testCase.Fixture.call();
            testCase.verifyTrue(r.finite);
            testCase.verifyTrue(r.symmetryResolution.swapTriggered);
            testCase.verifyFalse(r.symmetryResolution.swapAccepted);
            testCase.verifyEqual(r.callCounts.map,5);
            testCase.verifyEqual(r.selectedIndex,2);
            testCase.verifyEqual(r.srSwapException.identifier,"refine_position_map:PositionAtRisCenter");
        end
        function finiteNonconvergedStillSelected(testCase)
            testCase.Fixture.setMode('nonconverged');r=testCase.Fixture.call();
            testCase.verifyTrue(r.finite);
            testCase.verifyFalse(r.converged);
            testCase.verifyEqual(r.selectedIndex,1);
            testCase.verifyEqual(r.map.exitflag,0);
        end
        function programErrorPropagates(testCase)
            testCase.Fixture.setMode('program');
            testCase.verifyError(@()testCase.Fixture.call(),'s4Fixture:ProgramError');
        end
        function invalidConfigurationPropagates(testCase)
            testCase.verifyError(@()testCase.Fixture.callBadPrior(),'baseline2S4:InvalidPrior');
        end
    end
end
