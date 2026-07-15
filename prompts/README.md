# Prompts

Prompt resources used for AI-assisted GEP annotation.

GEPs were generated from the PedSTAR single-cell atlas using CSI-GEP, an NMF-based workflow for identifying robust gene expression programs.

Each GEP was summarized as a structured user prompt containing top loading genes, enrichment results, SingleR cell-type scores, NumBat malignancy probabilities, and disease-of-origin information for top-activating cells.

The AI agent analyzed each GEP independently using a shared system prompt and returned a label, justification, and confidence score.

After individual annotation, GEP labels and justifications were reviewed together in a second meta-annotation step to identify broader biological themes and experimental hypotheses.

- `system_prompt.txt`: model role, context, task, and output format.
- `user_prompt_example_gep1.txt`: example GEP-specific user prompt.
- `experimental_plan_by_ai_agent.txt`: prompt and output for downstream experimental planning.
