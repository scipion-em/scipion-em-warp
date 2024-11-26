from pyworkflow.wizard import Wizard
from warp.protocols import ProtWarpTSMotionCorr


class ProtWarpMaxResolutionWizard(Wizard):
    _targets = [(ProtWarpTSMotionCorr, ['range_max'])]

    @staticmethod
    def _getRangeHighResolution(protocol):
        range_high = 4

        if protocol.inputTSMovies.hasValue():
            sr = protocol.inputTSMovies.get().getSamplingRate()
            return round(sr * 2 + 0.1, 2)

        return range_high

    def show(self, form, *args):
        form.setVar('range_max', self._getRangeHighResolution(form.protocol))
