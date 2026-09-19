# AnHai2026-JCAM-WidthCrack

# Nhóm 2 — Width of Compact Sets & Ứng dụng đo độ rộng vết nứt

Repo triển khai lại các nội dung trong tài liệu:

> **An, P.T. & Hai, N.N. (2026).**"Widths of compact sets in normed spaces and an application in measuring crack widths." 
> Journal of Computational and Applied Mathematics **485**, 117498

---

## 1. Giới thiệu đề tài

Độ rộng của một hình chữ nhật là độ dài cạnh ngắn hơn — điều này hiển nhiên. Nhưng **độ rộng của một vết nứt trên tường, một ổ gà trên đường, hay một khối u thì được định nghĩa thế nào?**

Cho tới thời điểm bài báo được viết, các công trình đo độ rộng vết nứt (Tang et al., Wu et al.) đều dùng những quy ước đo mang tính kỹ thuật — theo pixel, theo hướng vuông góc trục vết nứt, theo bounding box — mà **không có định nghĩa toán học thống nhất**, nên hai phương pháp khác nhau cho hai kết quả khác nhau mà không có cơ sở so sánh. Ứng viên có sẵn trong hình học lồi là *mean width*, nhưng đại lượng này phụ thuộc chặt vào chu vi bao lồi — với một vết nứt dài, mảnh, ngoằn ngoèo thì bao lồi gần như không mang thông tin gì về độ rộng thực.

Bài báo đề xuất một định nghĩa **cục bộ** (tại từng điểm) thay vì toàn cục: độ rộng của $S$ tại $p$ là độ dài **dây cung ngắn nhất qua $p$** nối hai điểm biên. Cách này vừa cho kết quả đúng trực giác với hình chữ nhật, vừa đo được độ rộng *tại từng vị trí* dọc vết nứt bằng cách chọn $G$ là skeleton (trục giữa) của nó.

Repo này trả lời bốn câu hỏi:

1. Trình bày và **triển khai lại Algorithm 1** của bài báo. Thuật toán có hội tụ về nghiệm chính xác không?
2. Xây dựng **điều kiện đủ** để Algorithm 1 hội tụ về nghiệm chính xác.
3. Với $G \subset S$, xây dựng **thuật toán tính $W(S;G)$** trong $\mathbb{R}^2$.
4. **Ứng dụng thực tế**: đo độ rộng vết nứt tự phục hồi trên bộ dữ liệu Jakubowski & Tomczak (2025).

---

## 2. Cơ sở lý thuyết (tóm tắt)

Ký hiệu: $S \subset \mathbb{R}^2$ là tập compact, $\partial S$ là biên của $S$, $\text{int}(S)$ là phần trong của $S$, và $[x,y]$ là đoạn thẳng nối hai điểm $x$ và $y$.

**Definition 1–2.**

$$W(p) := \inf\{\|u-v\| : u,v \in \partial S,\ p \in [u,v]\}, \qquad W(G) := \sup_{p \in G} W(p)$$

Nếu $p \in \partial S$ thì $W(p)=0$; nếu $p \in \text{int}(S)$ thì $W(p) > 0$.

**Proposition 1.** Infimum trên **đạt được** — tồn tại $u^{\ast},v^{\ast} \in \partial S$ sao cho $\|u^{\ast}-v^{\ast}\| = W(p)$. 

**Theorem 1.** $x \mapsto W(x)$ **luôn nửa liên tục dưới (lsc)** trên $S$; liên tục tại $p \in \text{int}(S)$ nếu mọi tia từ $p$ cắt $\partial S$ tại đúng một điểm — tức không tồn tại hai điểm biên phân biệt $w, w' \in \partial S$ trên cùng một tia sao cho $w$ nằm giữa $p$ và $w'$, và đoạn thẳng nối $p$ với $w'$ nằm trọn trong $S$. 

**Example 3** trong bài báo cho một đa giác cụ thể mà tại đó $W$ **gián đoạn ở một điểm trong**, và tồn tại $G$ compact, $G \subset \text{int}(S)$ sao cho **không có điểm nào của $G$ đạt được $W(G)$** — tức supremum trong Definition 2 không phải lúc nào cũng đạt. 

**Theorem 2** (liên tục theo metric Hausdorff $d_{\mathcal H}$): $A \mapsto W(A)$ luôn lsc; nếu $K$ đóng, $K \subset \text{int}(S)$ và $W$ liên tục trên $K$ thì $A \mapsto W(A)$ liên tục đều trên các tập con đóng của $K$.

**Algorithm 1** (Section 5.2 của bài báo). Ý tưởng: $W(p)$ là cực tiểu theo *hướng* $d$ của độ dài dây cung qua $p$; thuật toán sinh $2^{n-1}$ hướng bằng cách chia đôi góc lặp, lấy đoạn ngắn nhất làm xấp xỉ $W_S^{(n)}(p)$.

**Proposition 7.** $\{W_S^{(n)}(p)\}$ là dãy giảm, hội tụ, và

$$\lim_{n\to\infty} W_S^{(n)}(p) \;\ge\; W_S(p)$$

— **chỉ chứng minh $\ge$, không chứng minh $=$**. Lý do: vì một đường thẳng qua $p$ theo hướng $\theta$ và hướng $\theta+\pi$ là cùng một đường (chord không phân biệt chiều), không gian hướng thực sự cần xét là $[0,\pi)$ chứ không phải $[0,2\pi)$. Tập hướng sinh bởi chia đôi góc là các góc nhị phân (dyadic) dạng $\theta_0 + k\pi/2^{n-1}$, trù mật trong $[0,\pi)$ nhưng nói chung **không chứa** hướng tối ưu $d^{\ast}$ tại mọi $n$ hữu hạn. Việc $\min$ trên dãy hướng trù mật có hội tụ về $\min$ toàn cục hay không phụ thuộc vào **tính liên tục** của hàm độ dài tại $d^{\ast}$.

**Trường hợp đa giác đơn — nghiệm chính xác.** Khi $S$ là đa giác, cặp $(u^*, v^*)$ đạt $W(p)$ nằm trên một cặp cạnh cụ thể; vì số cặp cạnh hữu hạn nên duyệt toàn bộ cặp và giải bài toán con "đoạn ngắn nhất qua $p$ trong một góc" (Problem 1, Appendix A — có nghiệm đóng khi $\theta=\pi/2$ hoặc khi $p$ nằm trên phân giác, nghiệm tổng quát qua phương trình bậc ba/Cardano) sẽ cho $W(p)$ **chính xác**.

---

## 3. Cấu trúc repo

| Thư mục | Nội dung | Câu hỏi | Chi tiết |
|---|---|---|---|
| `task1_algorithm1/` | Triển khai lại Algorithm 1, phân tích tính hội tụ | (1) | [EXPLANATION.md](task1_algorithm1/EXPLANATION.md) |
| `task2_convergence_conditions/` | Điều kiện đủ để Algorithm 1 hội tụ về nghiệm chính xác | (2) | [conditions.md](task2_convergence_conditions/conditions.md) |
| `task3_width_of_subset/` | Thuật toán tính $W(S;G)$ với $G \subset S$ trong $\mathbb{R}^2$ | (3) | [pseudocode.md](task3_width_of_subset/pseudocode.md) |
| `task4_selfhealing_application/` | Áp dụng đo độ rộng vết nứt tự phục hồi trên dataset Jakubowski & Tomczak (2025) | (4) | [README.md](task4_selfhealing_application/README.md) |
| `common/` | Hàm hình học & $W(p)$ dùng chung cho các task | — | — |
| `docs/` | Bài báo gốc, báo cáo cuối kỳ | — | — |

---

## 4. Cài đặt

```bash
python -m venv venv
source venv/bin/activate       # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

---

## 5. Chạy nhanh

Chạy các lệnh dưới đây từ thư mục gốc của repo (đã kích hoạt `venv`):

```bash
# (1) Demo Algorithm 1: hội tụ trên hình chữ nhật, so sánh với nghiệm chính xác
python task1_algorithm1/algorithm1.py

# (2) Kiểm chứng điều kiện hội tụ: S lồi (hội tụ đúng mọi θ0) vs. S không lồi (phụ thuộc θ0)
python task2_convergence_conditions/verify_conditions.py

# (3) Demo tính W(S;G) trên tập con G (đoạn thẳng qua Example 3)
python task3_width_of_subset/set_width_algorithm.py

# (4) Áp dụng lên ảnh vết nứt tự phục hồi
#     4a. Nếu chưa có dữ liệu thật, tạo ảnh mô phỏng để chạy thử pipeline:
python task4_selfhealing_application/make_synthetic_sample.py
#     4b. Chạy Algorithm 1 trên toàn bộ skeleton của từng ảnh:
python task4_selfhealing_application/run_algorithm1_on_cracks.py
#     4c. Phân tích xu hướng phục hồi theo thời gian + xuất biểu đồ:
python task4_selfhealing_application/healing_progress_analysis.py
```

**Kết quả kỳ vọng sau khi chạy:**

| Lệnh | File kết quả | Nội dung |
|---|---|---|
| `task1_algorithm1/algorithm1.py` | `task1_algorithm1/convergence_rectangle_center.png` | Đồ thị $W_S^{(n)}(p)$ hội tụ về $W(S)=4$ khi $n$ tăng |
| `task4_.../run_algorithm1_on_cracks.py` | `task4_selfhealing_application/results/width_summary.csv` | Bảng $W(S;G)$ theo từng ảnh/mốc thời gian |
| `task4_.../healing_progress_analysis.py` | `task4_selfhealing_application/results/figures/healing_curve.png` | Đường cong độ rộng vết nứt giảm dần theo thời gian |

---

## 6. Hướng tiếp cận từng task

**Task 1.** Biểu diễn $S$ bằng danh sách đỉnh đa giác, cài `ray_polygon_intersection` để tìm giao tia–biên, sinh $2^{n-1}$ hướng theo Algorithm 1, lấy min. Đối chiếu với nghiệm chính xác (duyệt cặp cạnh, mục 2) để đo sai số $W_S^{(n)}(p) - W_S(p)$ theo $n$ và kết luận về tính hội tụ. 

**Task 2.** Chứng minh: nếu hàm độ dài theo hướng liên tục tại hướng tối ưu $d^{\ast}$ thì $\lim_n W_S^{(n)}(p) = W_S(p)$. Giả thiết đủ để kiểm tra: $S$ lồi (Corollary 2 tự động thỏa điều kiện Theorem 1), hoặc $S$ *star-shaped* đối với $p$, hoặc yếu hơn — liên tục cục bộ quanh $d^{\ast}$. 

**Task 3.** Rời rạc hóa $G$ thành lưới điểm, tính $W(p_k)$ từng điểm (dùng Algorithm 1 hoặc nghiệm chính xác nếu $S$ là đa giác), lấy max theo Proposition 4. Dùng Theorem 2(b) để chặn sai số rời rạc hóa khi $G \subset \text{int}(S)$.

**Task 4.** Pipeline: ảnh → phân đoạn vết nứt ($S$) → trích skeleton ($G$) → $W(p)$ tại các điểm trên skeleton → $W(G)$. So sánh $W(G)$ giữa các mốc thời gian để định lượng tiến trình tự phục hồi. 

---

## 7. Dữ liệu tham khảo (Task 4)

>Jakubowski, J. & Tomczak, K. (2025). *Dataset for developing deep learning models to assess crack width and self-healing progress in concrete.* Scientific Data, **12**, 165.
https://www.researchgate.net/publication/388460885

>Bộ ảnh vết nứt để kiểm thử Algorithm 1 (từ bài báo gốc): https://drive.google.com/drive/folders/1zYe9Q0QQZGqB59Yl-fxw2EfIwgAr_a7n.
## 8. Tài liệu tham khảo chính

1. An, P.T., Hai, N.N. (2026). *Widths of compact sets in normed spaces and an application in measuring crack widths.* J. Comput. Appl. Math. 485, 117498.
2. An, P.T., Phuong, N.M., Quang, D.L., Trang, H.L. (2024). *A modified pothole detection approach to capture their width.* Phys. Scr. 99(4), 045007.
3. Melzak, Z.A. (1960). *A note on sets of constant width.* Proc. Am. Math. Soc. 11, 493–497.
4. Klette, R., Rosenfeld, A. (2004). *Digital Geometry.* Morgan Kaufmann.

---

## 9. Nhóm thực hiện

Nhóm 2 — Đồ án môn học, dựa trên tài liệu *AnHai2026-JCAM-WidthCrack.pdf*.

---
