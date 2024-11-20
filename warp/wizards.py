from pyworkflow.wizard import Wizard
from warp.protocols import ProtWarpTSMotionCorr


class ProtWarpMaxResolutionWizard(Wizard):
    _targets = [(ProtWarpTSMotionCorr, ['range_max'])]

    def _getRangeHighResolution(self, protocol):
        range_high = 4

        if protocol.inputSet.hasValue():
            sr = protocol.inputSet.get().getSamplingRate()
            return round(sr * 2 + 0.1, 2)

        return range_high

    def show(self, form, *args):
        form.setVar('range_high', self._getRangeHighResolution(form.protocol))
