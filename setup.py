from setuptools import setup, find_packages

setup(
    name="turbomole_analyzer",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["pydantic", "numpy", "pandas"],
    entry_points={
        "console_scripts": [
            # Updated the executable command name to match the script filename
            "turbomole_analyzer_cmd=turbomole_analyzer.turbomole_analyzer_cmd:main",
        ],
    },
)