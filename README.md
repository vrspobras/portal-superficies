# Portal de Superfícies — Codespaces

Portal local/remoto para SAU3 e SAU4, com leitura de TXT e geração de TIN em Python.

## Rodar no GitHub Codespaces

1. Crie um repositório no GitHub e envie estes arquivos.
2. Abra o repositório e use **Code → Codespaces → Create codespace on main**.
3. Aguarde o Codespace criar o ambiente e instalar `numpy` e `scipy`.
4. No terminal do Codespace, execute:

```bash
python app.py
```

5. O Codespaces deve encaminhar a porta 8000; abra o endereço mostrado na aba **Ports**.

## Teste de saúde

Abra `/api/test` e confirme `{"ok": true}`.

## Observação

O Codespace é uma máquina virtual remota criada para o projeto. O navegador do seu computador acessa o portal pelo encaminhamento da porta; o `python.exe` não é executado na sua máquina corporativa.
