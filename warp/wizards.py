from pyworkflow.wizard import Wizard
from warp.protocols import ProtWarpTSCtfEstimation


class ProtWarpMaxResolutionWizard(Wizard):
    _targets = [(ProtWarpTSCtfEstimation, ['range_high'])]

    def _getRangeHighResolution(self, protocol):
        range_high = 4

        if protocol.inputSet.hasValue():
            sr = protocol.inputSet.get().getSamplingRate()
            return sr * 2 + 0.1

        return range_high

    def show(self, form, *args):
        form.setVar('range_high', self._getRangeHighResolution(form.protocol))
