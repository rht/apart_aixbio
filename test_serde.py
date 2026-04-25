from langgraph.checkpoint.memory import MemorySaver
from langgraph.checkpoint.serde.jsonplus import JsonPlusSerializer

serde = JsonPlusSerializer(allowed_msgpack_modules=[
    ("aixbio.models.protein", "Chain"),
    ("aixbio.models.protein", "ProteinRecord"),
    ("aixbio.models.audit", "AgentDecision"),
])
m = MemorySaver(serde=serde)
print("MemorySaver initialized with serde successfully")
