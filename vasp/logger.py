# ###################################################################
# Logger for handling information, warning and debugging
# ###################################################################
import logging
import sys

logging.basicConfig(stream=sys.stdout,
                    format=('%(levelname)-8s [[%(pathname)s::%(lineno)d]'
                            '[%(funcName)s]]: %(message)s'),
                    level=logging.INFO)

log = logging.getLogger('Vasp')
