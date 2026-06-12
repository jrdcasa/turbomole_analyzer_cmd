from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

class BaseParser(ABC):
    @abstractmethod
    def parse(self, file_path: Path) -> Any:
        pass
